// Wright-Fisher simulation for a fwdpy11::DiploidPopulation with
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
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <fwdpp/diploid.hh>
#include <fwdpp/simparams.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>
#include <fwdpy11/evolvets/evolve_generation_ts.hpp>
#include <fwdpy11/evolvets/simplify_tables.hpp>
#include <fwdpy11/evolvets/sample_recorder_types.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpp/ts/recycling.hpp>
#include "util.hpp"
#include "diploid_pop_fitness.hpp"
#include "index_and_count_mutations.hpp"
#include "cleanup_metadata.hpp"
#include "track_mutation_counts.hpp"
#include "remove_extinct_mutations.hpp"
#include "track_ancestral_counts.hpp"
#include "remove_extinct_genomes.hpp"

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

void
apply_treseq_resetting_of_ancient_samples(
    const fwdpy11::DiploidPopulation_temporal_sampler &recorder,
    fwdpy11::DiploidPopulation &pop)
{
    recorder(pop);
    if (!pop.tables.preserved_nodes.empty())
        {
            pop.tables.preserved_nodes.clear();
            pop.ancient_sample_metadata.clear();
            pop.ancient_sample_genetic_value_matrix.clear();
        }
}

std::pair<std::vector<fwdpp::ts::TS_NODE_INT>, std::vector<std::size_t>>
simplification(
    bool preserve_selected_fixations, bool simulating_neutral_variants,
    bool suppress_edge_table_indexing,
    bool reset_treeseqs_to_alive_nodes_after_simplification,
    const fwdpy11::DiploidPopulation_temporal_sampler &post_simplification_recorder,
    fwdpp::ts::table_simplifier &simplifier, fwdpy11::DiploidPopulation &pop)
{
    auto simplification_rv = fwdpy11::simplify_tables(
        pop, pop.mcounts_from_preserved_nodes, pop.tables, simplifier,
        preserve_selected_fixations, simulating_neutral_variants,
        suppress_edge_table_indexing);
    if (pop.mcounts.size() != pop.mcounts_from_preserved_nodes.size())
        {
            throw std::runtime_error("evolvets: count vector size mismatch after "
                                     "simplification");
        }
    remap_metadata(pop.ancient_sample_metadata, simplification_rv.first);
    remap_metadata(pop.diploid_metadata, simplification_rv.first);
    if (reset_treeseqs_to_alive_nodes_after_simplification == true)
        {
            apply_treseq_resetting_of_ancient_samples(post_simplification_recorder, pop);
        }
    return simplification_rv;
}

void
final_population_cleanup(
    bool suppress_edge_table_indexing, bool preserve_selected_fixations,
    bool remove_extinct_mutations_at_finish, std::uint32_t last_preserved_generation,
    const std::vector<std::uint32_t> &last_preserved_generation_counts,
    fwdpy11::DiploidPopulation &pop)
{
    index_and_count_mutations(suppress_edge_table_indexing, pop);
    if (pop.generation == last_preserved_generation)
        {
            std::transform(begin(pop.mcounts_from_preserved_nodes),
                           end(pop.mcounts_from_preserved_nodes),
                           begin(last_preserved_generation_counts),
                           begin(pop.mcounts_from_preserved_nodes),
                           std::minus<std::uint32_t>());
        }
    if (!preserve_selected_fixations)
        {
            auto itr = std::remove_if(
                pop.tables.mutation_table.begin(), pop.tables.mutation_table.end(),
                [&pop](const fwdpp::ts::mutation_record &mr) {
                    return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                           && pop.mcounts_from_preserved_nodes[mr.key] == 0;
                });
            auto d = std::distance(itr, end(pop.tables.mutation_table));
            pop.tables.mutation_table.erase(itr, end(pop.tables.mutation_table));
            if (d)
                {
                    pop.tables.rebuild_site_table();
                }
            fwdpp::ts::remove_fixations_from_haploid_genomes(
                pop.haploid_genomes, pop.mutations, pop.mcounts,
                pop.mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                preserve_selected_fixations);
            // NOTE: this is hacky and should be better-handled upstream
            for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i] == 2 * pop.N
                        && pop.mcounts_from_preserved_nodes[i] == 0)
                        {
                            pop.mcounts[i] = 0;
                        }
                }
        }
    cleanup_metadata(pop.tables, pop.generation, pop.ancient_sample_metadata);
    if (remove_extinct_mutations_at_finish)
        {
            remove_extinct_mutations(pop);
        }
    remove_extinct_genomes(pop);
}

void
evolve_with_tree_sequences(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    fwdpy11::SampleRecorder &sr, const unsigned simplification_interval,
    ddemog::DiscreteDemography &demography, const std::uint32_t simlen,
    const double mu_neutral, const double mu_selected,
    const fwdpy11::MutationRegions &mmodel, const fwdpy11::GeneticMap &rmodel,
    // NOTE: gvalue_pointers is a change in 0.6.0,
    // and the object holds non-const bare pointers
    // to objects owned by Python.
    const fwdpy11::dgvalue_pointer_vector_ &gvalue_pointers,
    fwdpy11::DiploidPopulation_sample_recorder recorder,
    std::function<bool(const fwdpy11::DiploidPopulation &, const bool)>
        &stopping_criteron,
    const double selfing_rate,
    // NOTE: this is the complement of what a user will input, which is "prune_selected"
    const bool preserve_selected_fixations, const bool suppress_edge_table_indexing,
    bool record_genotype_matrix, const bool track_mutation_counts_during_sim,
    const bool remove_extinct_mutations_at_finish,
    const bool reset_treeseqs_to_alive_nodes_after_simplification,
    const fwdpy11::DiploidPopulation_temporal_sampler &post_simplification_recorder)
{
    if (gvalue_pointers.genetic_values.empty())
        {
            throw std::invalid_argument("empty list of genetic values");
        }
    //validate the input params
    if (pop.tables.genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument(
                "Population is not initialized with tree sequence support");
        }
    if (!std::isfinite(mu_selected))
        {
            throw std::invalid_argument("selected mutation rate is not finite");
        }
    if (mu_selected < 0.0)
        {
            throw std::invalid_argument("selected mutation rate must be non-negative");
        }
    if (!std::isfinite(mu_neutral))
        {
            throw std::invalid_argument("neutral mutation rate is not finite");
        }
    if (mu_neutral < 0.0)
        {
            throw std::invalid_argument("neutral mutation rate must be non-negative");
        }
    if (mu_neutral + mu_selected > 0.0 && mmodel.weights.empty())
        {
            throw std::invalid_argument(
                "nonzero mutation rate incompatible with empty regions");
        }
    if (!simlen)
        {
            throw std::invalid_argument("simulation length must be > 0");
        }
    if (pop.tables.node_table.empty())
        {
            throw std::invalid_argument("node table is not initialized");
        }

    // Set up discrete demography types. New in 0.6.0
    auto current_demographic_state = ddemog::initialize_model_state(
        pop.generation, pop.diploid_metadata, demography);
    if (gvalue_pointers.genetic_values.size()
        > static_cast<std::size_t>(current_demographic_state->maxdemes))
        {
            throw std::invalid_argument(
                "list of genetic values is longer than maxdemes");
        }
    if (static_cast<std::size_t>(current_demographic_state->maxdemes) > 1
        && gvalue_pointers.genetic_values.size() > 1
        && gvalue_pointers.genetic_values.size()
               < static_cast<std::size_t>(current_demographic_state->maxdemes))
        {
            throw std::invalid_argument("too few genetic value objects");
        }

    double total_mutation_rate = mu_neutral + mu_selected;
    const auto bound_mmodel = [&rng, &mmodel, &pop, total_mutation_rate](
                                  fwdpp::flagged_mutation_queue &recycling_bin,
                                  std::vector<fwdpy11::Mutation> &mutations) {
        std::vector<fwdpp::uint_t> rv;
        unsigned nmuts = gsl_ran_poisson(rng.get(), total_mutation_rate);
        for (unsigned i = 0; i < nmuts; ++i)
            {
                std::size_t x = gsl_ran_discrete(rng.get(), mmodel.lookup.get());
                auto key = mmodel.regions[x]->operator()(
                    recycling_bin, mutations, pop.mut_lookup, pop.generation, rng);
                rv.push_back(key);
            }
        std::sort(begin(rv), end(rv),
                  [&mutations](const fwdpp::uint_t a, const fwdpp::uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    };

    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng); };

    auto genetics = fwdpp::make_genetic_parameters(gvalue_pointers.genetic_values,
                                                   std::move(bound_mmodel),
                                                   std::move(bound_rmodel));
    std::vector<std::size_t> deme_to_gvalue_map(current_demographic_state->maxdemes, 0);
    if (genetics.gvalue.size() > 1)
        {
            std::iota(begin(deme_to_gvalue_map), end(deme_to_gvalue_map), 0);
        }
    // A stateful fitness model will need its data up-to-date,
    // so we must call update(...) prior to calculating fitness,
    // else bad stuff like segfaults could happen.
    for (auto &i : genetics.gvalue)
        {
            i->update(pop);
        }
    std::vector<fwdpy11::DiploidMetadata> offspring_metadata(pop.diploid_metadata);
    std::vector<fwdpy11::DiploidGenotype> offspring;
    std::vector<double> new_diploid_gvalues;
    calculate_diploid_fitness(rng, pop, genetics.gvalue, deme_to_gvalue_map,
                              offspring_metadata, new_diploid_gvalues,
                              record_genotype_matrix);
    pop.diploid_metadata.swap(offspring_metadata);

    if (!pop.mutations.empty())
        {
            // It is possible that pop already has a tree sequence
            // containing neutral variants not in haploid_genome objects.
            // To correctly handle this case, we build the recycling bin
            // from any elements in pop.mutations not in the mutation table.
            if (!pop.tables.node_table.empty())
                {
                    std::vector<std::size_t> indexes(pop.mutations.size());
                    std::iota(begin(indexes), end(indexes), 0);
                    for (auto &m : pop.tables.mutation_table)
                        {
                            indexes[m.key] = std::numeric_limits<std::size_t>::max();
                        }
                    std::queue<std::size_t> can_recycle;
                    for (auto i : indexes)
                        {
                            if (i != std::numeric_limits<std::size_t>::max())
                                {
                                    can_recycle.push(i);
                                }
                        }
                    genetics.mutation_recycling_bin
                        = fwdpp::flagged_mutation_queue(std::move(can_recycle));
                }
            else // Assume there are no tree sequence data
                {
                    // Then we assume pop exists in an "already simulated"
                    // state and is properly-book-kept
                    genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                        pop.mcounts, pop.mcounts_from_preserved_nodes);
                }
        }

    fwdpp::ts::TS_NODE_INT next_index = pop.tables.node_table.size();
    bool simplified = false;
    fwdpp::ts::table_simplifier simplifier(pop.tables.genome_length());
    bool stopping_criteron_met = false;
    const bool simulating_neutral_variants = (mu_neutral > 0.0) ? true : false;
    std::pair<std::vector<fwdpp::ts::TS_NODE_INT>, std::vector<std::size_t>>
        simplification_rv;
    std::uint32_t last_preserved_generation = std::numeric_limits<std::uint32_t>::max();
    decltype(pop.mcounts) last_preserved_generation_counts;
    for (std::uint32_t gen = 0; gen < simlen && !stopping_criteron_met; ++gen)
        {
            // NOTE: demographic changes are applied with respect
            // to a PARENTAL generation, so we call this...
            ddemog::update_demography_manager(rng, pop.generation, pop.diploid_metadata,
                                              demography, current_demographic_state);

            if (current_demographic_state->will_go_globally_extinct() == true)
                {
                    simplification_rv = simplification(
                        preserve_selected_fixations, simulating_neutral_variants,
                        suppress_edge_table_indexing,
                        reset_treeseqs_to_alive_nodes_after_simplification,
                        post_simplification_recorder, simplifier, pop);
                    final_population_cleanup(
                        suppress_edge_table_indexing, preserve_selected_fixations,
                        remove_extinct_mutations_at_finish, last_preserved_generation,
                        last_preserved_generation_counts, pop);
                    std::ostringstream o;
                    o << "extinction at time " << pop.generation;
                    throw ddemog::GlobalExtinction(o.str());
                }
            // ...before updating this:
            ++pop.generation;
            fwdpy11::evolve_generation_ts(rng, pop, genetics, current_demographic_state,
                                          pop.generation, pop.tables, offspring,
                                          offspring_metadata, next_index);
            // TODO: abstract out these steps into a "cleanup_pop" function
            // TODO: this needs cleanup in 0.7.0, so that parental genotypes
            // can affect things.  Right now, we are requiring the offspring
            // genotypes to be in pop.diploids in the gvalue calculations.
            pop.diploids.swap(offspring);
            // TODO: deal with random effects
            for (auto &i : genetics.gvalue)
                {
                    i->update(pop);
                }

            calculate_diploid_fitness(rng, pop, genetics.gvalue, deme_to_gvalue_map,
                                      offspring_metadata, new_diploid_gvalues,
                                      record_genotype_matrix);
            // NOTE: this may need revisiting in 0.7.0
            pop.genetic_value_matrix.swap(new_diploid_gvalues);
            // TODO: abstract out these steps into a "cleanup_pop" function
            pop.diploid_metadata.swap(offspring_metadata);
            pop.N = static_cast<std::uint32_t>(pop.diploid_metadata.size());

            if (gen > 0 && gen % simplification_interval == 0.0)
                {
                    simplification_rv = simplification(
                        preserve_selected_fixations, simulating_neutral_variants,
                        suppress_edge_table_indexing,
                        reset_treeseqs_to_alive_nodes_after_simplification,
                        post_simplification_recorder, simplifier, pop);
                    simplified = true;
                }
            else
                {
                    simplified = false;
                }
            if (pop.tables.num_nodes()
                >= std::numeric_limits<fwdpp::ts::TS_NODE_INT>::max() - 1)
                {
                    throw std::runtime_error("range error for node labels");
                }
            next_index = pop.tables.num_nodes();
            if (track_mutation_counts_during_sim)
                {
                    track_mutation_counts(pop, simplified, suppress_edge_table_indexing);
                }
            // The user may now analyze the pop'n and record ancient samples
            recorder(pop, sr);
            if (simplified)
                {
                    if (suppress_edge_table_indexing == false)
                        {
                            // Behavior change in 0.5.3: set all fixation counts to 0
                            // to flag for recycling if possible.
                            // NOTE: this may slow things down a touch?
                            if (preserve_selected_fixations == false)
                                {
                                    // b/c neutral mutations not in genomes!
                                    fwdpp::ts::remove_fixations_from_haploid_genomes(
                                        pop.haploid_genomes, pop.mutations, pop.mcounts,
                                        pop.mcounts_from_preserved_nodes,
                                        2 * pop.diploids.size(),
                                        preserve_selected_fixations);
                                }
                            for (auto &i : simplification_rv.second)
                                {
                                    if (pop.mcounts[i] == 2 * pop.diploids.size()
                                        && pop.mcounts_from_preserved_nodes[i] == 0)
                                        {
                                            if (pop.mutations[i].neutral
                                                || !preserve_selected_fixations)
                                                {
                                                    // flag variant for recycling
                                                    pop.mcounts[i] = 0;
                                                    // flag item for removal from return value,
                                                    // as mutation is no longer considered "preserved"
                                                    i = std::numeric_limits<
                                                        std::size_t>::max();
                                                }
                                        }
                                }
                            simplification_rv.second.erase(
                                std::remove(begin(simplification_rv.second),
                                            end(simplification_rv.second),
                                            std::numeric_limits<std::size_t>::max()),
                                end(simplification_rv.second));
                            genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                pop.mcounts, pop.mcounts_from_preserved_nodes);
                        }
                    else
                        {
                            genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                simplification_rv.second, pop.mutations.size());
                        }
                }

            // TODO: deal with the result of the recorder populating sr
            if (!sr.samples.empty())
                {
                    for (auto i : sr.samples)
                        {
                            if (i >= pop.N)
                                {
                                    throw std::invalid_argument(
                                        "ancient sample index greater than "
                                        "current population size");
                                }
                            auto x = fwdpy11::parent_nodes_from_metadata(
                                i, pop.diploid_metadata, 0);
                            pop.tables.preserved_nodes.push_back(x.first);
                            pop.tables.preserved_nodes.push_back(x.second);

                            // Record the metadata for this individual
                            pop.ancient_sample_metadata.push_back(
                                pop.diploid_metadata[i]);
                            // Update the genotype matrix w.r.to
                            // the new ancient samples
                            if (!pop.genetic_value_matrix.empty())
                                {
                                    auto offset = i * genetics.gvalue[0]->total_dim;
                                    pop.ancient_sample_genetic_value_matrix.insert(
                                        end(pop.ancient_sample_genetic_value_matrix),
                                        begin(pop.genetic_value_matrix) + offset,
                                        begin(pop.genetic_value_matrix) + offset
                                            + genetics.gvalue[0]->total_dim);
                                }
                        }
                    track_ancestral_counts(sr.samples, &last_preserved_generation,
                                           last_preserved_generation_counts, pop);
                    // Finally, clear the input
                    sr.samples.clear();
                }
            stopping_criteron_met = stopping_criteron(pop, simplified);
        }

    // NOTE: if tables.preserved_nodes overlaps with samples,
    // then simplification throws an error. But, since it is annoying
    // for a user to have to remember not to do that, we filter the list
    // here
    pop.fill_alive_nodes();
    std::sort(begin(pop.alive_nodes), end(pop.alive_nodes));
    auto itr = std::remove_if(
        begin(pop.tables.preserved_nodes), end(pop.tables.preserved_nodes),
        [&pop](const fwdpp::ts::TS_NODE_INT l) {
            return std::binary_search(begin(pop.alive_nodes), end(pop.alive_nodes), l);
        });
    pop.tables.preserved_nodes.erase(itr, end(pop.tables.preserved_nodes));
    pop.alive_nodes.clear();

    if (!simplified)
        {
            simplification_rv = simplification(
                preserve_selected_fixations, simulating_neutral_variants,
                suppress_edge_table_indexing,
                reset_treeseqs_to_alive_nodes_after_simplification,
                post_simplification_recorder, simplifier, pop);
        }
    final_population_cleanup(suppress_edge_table_indexing, preserve_selected_fixations,
                             remove_extinct_mutations_at_finish,
                             last_preserved_generation, last_preserved_generation_counts,
                             pop);
    ddemog::save_model_state(current_demographic_state, demography);
}

void
init_evolve_with_tree_sequences(py::module &m)
{
    m.def("evolve_with_tree_sequences", &evolve_with_tree_sequences);
}
