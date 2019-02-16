// Wright-Fisher simulation for a fwdpy11::SlocusPop with
// tree sequences.
//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <queue>
#include <cmath>
#include <stdexcept>
#include <fwdpp/diploid.hh>
#include <fwdpp/simparams.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/samplers.hpp>
//#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/evolvets/evolve_generation_ts.hpp>
#include <fwdpy11/evolvets/simplify_tables.hpp>
#include <fwdpy11/evolvets/sample_recorder_types.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>

namespace py = pybind11;

inline void
resize_genotype_matrix(std::vector<double> &new_diploid_gvalues,
                       std::size_t newsize, std::true_type)
{
    new_diploid_gvalues.resize(newsize);
}

inline void
resize_genotype_matrix(std::vector<double> & /*new_diploid_gvalues*/,
                       std::size_t /*newsize*/, std::false_type)
{
}

template <typename gvalue_fxn>
inline void
copy_genetic_values(double *beg, const gvalue_fxn &genetic_value_fxn,
                    std::true_type)
{
    std::copy(begin(genetic_value_fxn.gvalues), end(genetic_value_fxn.gvalues),
              beg);
}

template <typename gvalue_fxn>
inline void
copy_genetic_values(double * /*beg*/, const gvalue_fxn & /*genetic_value_fxn*/,
                    std::false_type)
{
}

template <typename update_genotype_matrix>
fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitness_details(
    const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
    const fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
    std::vector<fwdpy11::DiploidMetadata> &new_metadata,
    std::vector<double> &new_diploid_gvalues, const update_genotype_matrix um)
{
    // Calculate parental fitnesses
    std::vector<double> parental_fitnesses(pop.diploids.size());
    double sum_parental_fitnesses = 0.0;
    new_metadata.resize(pop.N);
    resize_genotype_matrix(new_diploid_gvalues,
                           pop.N * genetic_value_fxn.total_dim, um);
    auto gvoffset = new_diploid_gvalues.data();
    for (std::size_t i = 0; i < pop.diploids.size();
         ++i, gvoffset += genetic_value_fxn.total_dim)
        {
            new_metadata[i] = pop.diploid_metadata[i];
            genetic_value_fxn(rng, i, pop, new_metadata[i]);
            copy_genetic_values(gvoffset, genetic_value_fxn, um);
            parental_fitnesses[i] = new_metadata[i].w;
            sum_parental_fitnesses += parental_fitnesses[i];
        }
    pop.diploid_metadata.swap(new_metadata);
    pop.genetic_value_matrix.swap(new_diploid_gvalues);
    // If the sum of parental fitnesses is not finite,
    // then the genetic value calculator returned a non-finite value/
    // Unfortunately, gsl_ran_discrete_preproc allows such values through
    // without raising an error, so we have to check things here.
    if (!std::isfinite(sum_parental_fitnesses))
        {
            throw std::runtime_error("non-finite fitnesses encountered");
        }

    auto rv = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(parental_fitnesses.size(),
                                 parental_fitnesses.data()));
    if (rv == nullptr)
        {
            // This is due to negative fitnesses
            throw std::runtime_error(
                "fitness lookup table could not be generated");
        }
    return rv;
}

std::function<fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
    const fwdpy11::GSLrng_t &g, fwdpy11::SlocusPop &,
    const fwdpy11::SlocusPopGeneticValue &,
    std::vector<fwdpy11::DiploidMetadata> &, std::vector<double> &)>
wrap_calculate_fitness(bool update_genotype_matrix)
{
    if (update_genotype_matrix)
        {
            return [](const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
                      const fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
                      std::vector<fwdpy11::DiploidMetadata> &new_metadata,
                      std::vector<double> &new_diploid_gvalues) {
                return calculate_fitness_details(
                    rng, pop, genetic_value_fxn, new_metadata,
                    new_diploid_gvalues, std::true_type());
            };
        }
    return [](const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
              const fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
              std::vector<fwdpy11::DiploidMetadata> &new_metadata,
              std::vector<double> &new_diploid_gvalues) {
        return calculate_fitness_details(rng, pop, genetic_value_fxn,
                                         new_metadata, new_diploid_gvalues,
                                         std::false_type());
    };
}

// TODO: put in header for reuse
template <typename poptype>
void
remap_ancient_samples(poptype &pop,
                      const std::vector<fwdpp::ts::TS_NODE_INT> &idmap)
{
    for (auto &a : pop.ancient_sample_records)
        {
            a.n1 = idmap[a.n1];
            a.n2 = idmap[a.n2];
            if (a.n1 == fwdpp::ts::TS_NULL_NODE
                || a.n2 == fwdpp::ts::TS_NULL_NODE)
                {
                    throw std::runtime_error(
                        "error simplifying with respect to ancient samples");
                }
        }
}

// TODO: allow for neutral mutations in the future
void
wfSlocusPop_ts(
    const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
    fwdpy11::samplerecorder &sr, const unsigned simplification_interval,
    py::array_t<std::uint32_t> popsizes, //const double mu_neutral,
    const double mu_selected, const double recrate,
    const fwdpy11::MutationRegions &mmodel,
    const fwdpy11::RecombinationRegions &rmodel,
    fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
    fwdpy11::SlocusPop_sample_recorder recorder, const double selfing_rate,
    // NOTE: this is the complement of what a user will input, which is "prune_selected"
    const bool preserve_selected_fixations,
    const bool suppress_edge_table_indexing, bool record_genotype_matrix)
{
    //validate the input params
    if (pop.tables.genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument(
                "Population is not initialized with tree sequence support");
        }
    if (!std::isfinite(mu_selected))
        {
            throw std::invalid_argument(
                "selected mutation rate is not finite");
        }
    if (mu_selected < 0.0)
        {
            throw std::invalid_argument(
                "selected mutation rate must be non-negative");
        }
    const std::uint32_t num_generations
        = static_cast<std::uint32_t>(popsizes.size());
    if (!num_generations)
        {
            throw std::invalid_argument("empty list of population sizes");
        }
    if (pop.tables.node_table.empty())
        {
            throw std::invalid_argument("node table is not initialized");
        }

    const auto bound_mmodel = [&rng, &mmodel, &pop, mu_selected](
                                  fwdpp::flagged_mutation_queue &recycling_bin,
                                  std::vector<fwdpy11::Mutation> &mutations) {
        std::vector<fwdpp::uint_t> rv;
        unsigned nmuts = gsl_ran_poisson(rng.get(), mu_selected);
        for (unsigned i = 0; i < nmuts; ++i)
            {
                std::size_t x
                    = gsl_ran_discrete(rng.get(), mmodel.lookup.get());
                auto key = mmodel.regions[x]->operator()(
                    recycling_bin, mutations, pop.mut_lookup, pop.generation,
                    rng);
                rv.push_back(key);
            }
        std::sort(begin(rv), end(rv),
                  [&mutations](const fwdpp::uint_t a, const fwdpp::uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    };

    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng); };

    auto genetics = fwdpp::make_genetic_parameters(
        &genetic_value_fxn, std::move(bound_mmodel),
        std::move(bound_rmodel));
    // A stateful fitness model will need its data up-to-date,
    // so we must call update(...) prior to calculating fitness,
    // else bad stuff like segfaults could happen.
    genetic_value_fxn.update(pop);
    std::vector<fwdpy11::DiploidMetadata> new_metadata(pop.N);
    std::vector<double> new_diploid_gvalues;
    auto calculate_fitness = wrap_calculate_fitness(record_genotype_matrix);
    auto lookup = calculate_fitness(rng, pop, genetic_value_fxn,
                                    new_metadata, new_diploid_gvalues);

    // Generate our fxns for picking parents

    // Because lambdas that capture by reference do a "late" binding of
    // params, this is safe w.r.to updating lookup after each generation.
    const auto pick_first_parent = [&rng, &lookup]() {
        return gsl_ran_discrete(rng.get(), lookup.get());
    };

    const auto pick_second_parent
        = [&rng, &lookup, selfing_rate](const std::size_t p1) {
              if (selfing_rate == 1.0
                  || (selfing_rate > 0.0
                      && gsl_rng_uniform(rng.get()) < selfing_rate))
                  {
                      return p1;
                  }
              return gsl_ran_discrete(rng.get(), lookup.get());
          };
    const auto generate_offspring_metadata
        = [&rng](fwdpy11::DiploidMetadata &offspring_metadata,
                 const std::size_t p1, const std::size_t p2,
                 const std::vector<fwdpy11::DiploidMetadata>
                     & /*parental_metadata*/) {
              offspring_metadata.deme = 0;
              offspring_metadata.parents[0] = p1;
              offspring_metadata.parents[1] = p2;
          };
    fwdpp::flagged_mutation_queue mutation_recycling_bin
        = fwdpp::empty_mutation_queue();
    fwdpp::ts::TS_NODE_INT first_parental_index = 0,
                           next_index = pop.tables.node_table.size();
    bool simplified = false;
    fwdpp::ts::table_simplifier simplifier(pop.tables.genome_length());
    std::vector<fwdpp::ts::TS_NODE_INT> ancient_samples;
    for (std::uint32_t gen = 0; gen < num_generations; ++gen)
        {
            ++pop.generation;
            const auto N_next = popsizes.at(gen);
            // TODO: can simplify function further b/c we are referring
            // to data that fwdpy11 Populations contain.
            fwdpy11::evolve_generation_ts(
                rng, pop, genetics, N_next, pick_first_parent,
                pick_second_parent, generate_offspring_metadata,
                pop.generation, pop.tables, first_parental_index, next_index);

            //N_next, mu_selected, pick_first_parent,
            //pick_second_parent, generate_offspring_metadata, bound_mmodel,
            //mutation_recycling_bin, bound_rmodel, pop.generation,
            //first_parental_index, next_index);

            pop.N = N_next;
            // TODO: deal with random effects
            genetic_value_fxn.update(pop);
            lookup = calculate_fitness(rng, pop, genetic_value_fxn,
                                       new_metadata, new_diploid_gvalues);
            if (gen > 0 && gen % simplification_interval == 0.0)
                {
                    // TODO: update this to allow neutral mutations to be simulated
                    auto rv = fwdpy11::simplify_tables(
                        pop, pop.mcounts_from_preserved_nodes, pop.tables,
                        simplifier, pop.tables.num_nodes() - 2 * pop.N,
                        2 * pop.N, preserve_selected_fixations, false,
                        suppress_edge_table_indexing);
                    if (suppress_edge_table_indexing == false)
                        {
                            mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                pop.mcounts, pop.mcounts_from_preserved_nodes);
                        }
                    else
                        {
                            mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                rv.second, pop.mutations.size());
                        }
                    simplified = true;
                    next_index = pop.tables.num_nodes();
                    first_parental_index = 0;
                    remap_ancient_samples(pop, rv.first);
                }
            else
                {
                    simplified = false;
                    first_parental_index = next_index;
                    next_index += 2 * pop.N;
                }
            // The user may now analyze the pop'n and record ancient samples
            recorder(pop, sr);
            // TODO: deal with the result of the recorder populating sr
            if (!sr.samples.empty())
                {
                    ancient_samples.clear();
                    for (auto i : sr.samples)
                        {
                            if (i >= pop.N)
                                {
                                    throw std::invalid_argument(
                                        "ancient sample index greater than "
                                        "current population size");
                                }
                            // Get the nodes
                            auto x = fwdpp::ts::get_parent_ids(
                                first_parental_index, i, 0);
                            ancient_samples.push_back(x.first);
                            ancient_samples.push_back(x.second);

                            // Record the metadata for this individual
                            pop.ancient_sample_metadata.push_back(
                                pop.diploid_metadata[i]);
                            // Update the genotype matrix w.r.to
                            // the new ancient samples
                            if (!pop.genetic_value_matrix.empty())
                                {
                                    pop.ancient_sample_genetic_value_matrix.insert(
                                        end(pop.ancient_sample_genetic_value_matrix),
                                        begin(pop.genetic_value_matrix) + i,
                                        begin(pop.genetic_value_matrix) + i
                                            + genetic_value_fxn.total_dim);
                                }
                            // Record the time and nodes for this individual
                            pop.ancient_sample_records.emplace_back(
                                fwdpy11::ancient_sample_record{
                                    static_cast<double>(pop.generation),
                                    x.first, x.second });
                        }
                    // NOTE: this can throw an exception
                    pop.tables.record_preserved_nodes(ancient_samples);
                    // Finally, clear the input
                    sr.samples.clear();
                }
        }
    if (!simplified)
        {
            // NOTE: if tables.preserved_nodes overlaps with samples,
            // then simplification throws an error. But, since it is annoying
            // for a user to have to remember not to do that, we filter the list
            // here
            if (std::any_of(pop.tables.preserved_nodes.rbegin(),
                            pop.tables.preserved_nodes.rend(),
                            [&pop](const fwdpp::ts::TS_NODE_INT l) {
                                return l >= pop.tables.num_nodes() - 2 * pop.N;
                            }))
                {
                    auto itr = std::remove_if(
                        begin(pop.tables.preserved_nodes),
                        end(pop.tables.preserved_nodes),
                        [&pop](const fwdpp::ts::TS_NODE_INT l) {
                            return l >= pop.tables.num_nodes() - 2 * pop.N;
                        });
                    pop.tables.preserved_nodes.erase(
                        itr, end(pop.tables.preserved_nodes));
                }

            // TODO: update this to allow neutral mutations to be simulated
            auto rv = fwdpy11::simplify_tables(
                pop, pop.mcounts_from_preserved_nodes, pop.tables, simplifier,
                pop.tables.num_nodes() - 2 * pop.N, 2 * pop.N,
                preserve_selected_fixations, false,
                suppress_edge_table_indexing);

            remap_ancient_samples(pop, rv.first);
        }
    if (suppress_edge_table_indexing == true)
        {
            pop.tables.build_indexes();
            std::vector<std::int32_t> samples(2 * pop.N);
            std::iota(samples.begin(), samples.end(), 0);
            fwdpp::ts::count_mutations(pop.tables, pop.mutations, samples,
                                       pop.mcounts,
                                       pop.mcounts_from_preserved_nodes);
        }
}

PYBIND11_MODULE(wright_fisher_slocus_ts, m)
{
    m.doc() = "Evolution under a Wright-Fisher model using tree sequences.";

    m.def("WFSlocusPop_ts", &wfSlocusPop_ts);
}
