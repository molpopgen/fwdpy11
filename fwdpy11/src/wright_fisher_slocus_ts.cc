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
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/evolvets/evolve_generation_ts.hpp>
#include <fwdpy11/evolvets/sample_recorder_types.hpp>

namespace py = pybind11;

fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitness(const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
                  const fwdpy11::SlocusPopGeneticValue &genetic_value_fxn)
{
    // Calculate parental fitnesses
    std::vector<double> parental_fitnesses(pop.diploids.size());
    double sum_parental_fitnesses = 0.0;
    for (std::size_t i = 0; i < pop.diploids.size(); ++i)
        {
            auto g = genetic_value_fxn(i, pop);
            pop.diploid_metadata[i].g = g;
            pop.diploid_metadata[i].e = genetic_value_fxn.noise(
                rng, pop.diploid_metadata[i],
                pop.diploid_metadata[i].parents[0],
                pop.diploid_metadata[i].parents[1], pop);
            pop.diploid_metadata[i].w
                = genetic_value_fxn.genetic_value_to_fitness(
                    pop.diploid_metadata[i]);
            parental_fitnesses[i] = pop.diploid_metadata[i].w;
            sum_parental_fitnesses += parental_fitnesses[i];
        }
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


// TODO: allow for neutral mutations in the future
void
wfSlocusPop_ts(
    const fwdpy11::GSLrng_t &rng, fwdpy11::SlocusPop &pop,
    py::array_t<std::uint32_t> popsizes, //const double mu_neutral,
    const double mu_selected, const double recrate,
    const fwdpp::extensions::discrete_mut_model<fwdpy11::SlocusPop::mcont_t>
        &mmodel,
    const fwdpp::extensions::discrete_rec_model &rmodel,
    fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
    fwdpy11::SlocusPop_temporal_sampler recorder, const double selfing_rate,
    const bool remove_selected_fixations)
{
    //validate the input params
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
    const auto bound_mmodel = fwdpp::extensions::bind_dmm(rng.get(), mmodel);
    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng.get()); };

    // A stateful fitness model will need its data up-to-date,
    // so we must call update(...) prior to calculating fitness,
    // else bad stuff like segfaults could happen.
    genetic_value_fxn.update(pop);
    auto lookup = calculate_fitness(rng, pop, genetic_value_fxn);

    // Generate our fxns for picking parents

    // Because lambdas that capture by reference do a "late" binding of
    // params, this is safe w.r.to updating lookup after each generation.
    const auto pick_first_parent = [&rng, &lookup]() {
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
}

PYBIND11_MODULE(wright_fisher_slocus_ts, m)
{
    m.doc() = "Evolution under a Wright-Fisher model using tree sequences.";

    m.def("WFSlocusPop_ts", &wfSlocusPop_ts);
}
