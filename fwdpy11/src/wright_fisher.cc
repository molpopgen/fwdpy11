// Wright-Fisher simulation for a fwdpy11::DiploidPopulation
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
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/evolve/DiploidPopulation_generation.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>

namespace py = pybind11;

fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitness(const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
                  const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
                  std::vector<fwdpy11::DiploidMetadata> &new_metadata)
{
    // Calculate parental fitnesses
    std::vector<double> parental_fitnesses(pop.diploids.size());
    double sum_parental_fitnesses = 0.0;
    new_metadata.resize(pop.N);
    for (std::size_t i = 0; i < pop.diploids.size(); ++i)
        {
            new_metadata[i] = pop.diploid_metadata[i];
            genetic_value_fxn(rng, i, pop, new_metadata[i]);
            parental_fitnesses[i] = new_metadata[i].w;
            sum_parental_fitnesses += parental_fitnesses[i];
        }
    pop.diploid_metadata.swap(new_metadata);
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

void
handle_fixations(const bool remove_selected_fixations,
                 const std::uint32_t N_next, fwdpy11::DiploidPopulation &pop)
{
    if (remove_selected_fixations)
        {
            fwdpp::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                                  pop.mcounts, 2 * N_next,
                                                  std::true_type());
        }
    else
        {
            fwdpp::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                                  pop.mcounts, 2 * N_next,
                                                  fwdpp::remove_neutral());
        }
    fwdpy11::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                              pop.mut_lookup, pop.mcounts, pop.generation,
                              2 * pop.N, remove_selected_fixations);
}

void
wfDiploidPopulation(const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
            py::array_t<std::uint32_t> popsizes, const double mu_neutral,
            const double mu_selected, const double recrate,
            const fwdpy11::MutationRegions &mmodel,
            const fwdpy11::GeneticMap &rmodel,
            fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
            fwdpy11::DiploidPopulation_temporal_sampler recorder,
            const double selfing_rate, const bool remove_selected_fixations)
{
    //validate the input params
    if (!std::isfinite(mu_neutral))
        {
            throw std::invalid_argument("neutral mutation rate is not finite");
        }
    if (!std::isfinite(mu_selected))
        {
            throw std::invalid_argument(
                "selected mutation rate is not finite");
        }
    if (mu_neutral < 0.0)
        {
            throw std::invalid_argument(
                "neutral mutation rate must be non-negative");
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

    // E[S_{2N}] I got the expression from Ewens.
    pop.mutations.reserve(std::ceil(
        std::log(2 * pop.N) * (4. * double(pop.N) * (mu_neutral + mu_selected))
        + 0.667 * (4. * double(pop.N) * (mu_neutral + mu_selected))));

    const auto bound_mmodel
        = [&rng, &mmodel, &pop](fwdpp::flagged_mutation_queue &recycling_bin,
                                std::vector<fwdpy11::Mutation> &mutations) {
              std::size_t x = gsl_ran_discrete(rng.get(), mmodel.lookup.get());
              return mmodel.regions[x]->operator()(recycling_bin, mutations,
                                                   pop.mut_lookup,
                                                   pop.generation, rng);
          };
    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng); };

    // A stateful fitness model will need its data up-to-date,
    // so we must call update(...) prior to calculating fitness,
    // else bad stuff like segfaults could happen.
    genetic_value_fxn.update(pop);
    std::vector<fwdpy11::DiploidMetadata> new_metadata(pop.N);
    auto lookup = calculate_fitness(rng, pop, genetic_value_fxn, new_metadata);

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
        = [](fwdpy11::DiploidMetadata &offspring_metadata,
             const std::size_t p1, const std::size_t p2,
             const std::vector<fwdpy11::DiploidMetadata>
                 & /*parental_metadata*/) {
              offspring_metadata.deme = 0;
              offspring_metadata.parents[0] = p1;
              offspring_metadata.parents[1] = p2;
              offspring_metadata.nodes[0] = offspring_metadata.nodes[1] = -1;
          };

    for (std::uint32_t gen = 0; gen < num_generations; ++gen)
        {
            ++pop.generation;
            const auto N_next = popsizes.at(gen);
            fwdpy11::evolve_generation(
                rng, pop, N_next, mu_neutral + mu_selected, bound_mmodel,
                bound_rmodel, pick_first_parent, pick_second_parent,
                generate_offspring_metadata);
            handle_fixations(remove_selected_fixations, N_next, pop);

            pop.N = N_next;
            // TODO: deal with random effects
            genetic_value_fxn.update(pop);
            lookup
                = calculate_fitness(rng, pop, genetic_value_fxn, new_metadata);
            recorder(pop); // The user may now analyze the pop'n
        }
}

PYBIND11_MODULE(_wright_fisher, m)
{
    m.doc() = "Evolution under a Wright-Fisher model.";

    m.def("WFDiploidPopulation", &wfDiploidPopulation);
}
