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
#include <fwdpy11/types.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/rules/wf_rules.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/evolve/slocuspop.hpp>
namespace py = pybind11;

template <typename bound_mmodels, typename bound_recmodels,
          typename mut_removal_policy, typename update_mut>
void
evolve_common(const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
              fwdpy11::wf_rules& rules, py::array_t<std::uint32_t> popsizes,
              const double mu_neutral, const double mu_selected,
              const bound_mmodels& mmodels, const bound_recmodels& recmap,
              fwdpy11::single_locus_fitness& fitness,
              fwdpy11::singlepop_temporal_sampler& recorder,
              const double selfing_rate, const mut_removal_policy& mp,
              const update_mut& um)
{
    auto generations = popsizes.size();
    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    auto wbar = rules.w(pop, fitness_callback);
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            const auto N_next = popsizes.at(generation);
            fwdpy11::evolve_generation(
                rng, pop, N_next, mu_neutral + mu_selected, mmodels,
                recmap,
                std::bind(&fwdpy11::wf_rules::pick1, &rules,
                          std::placeholders::_1, std::placeholders::_2),
                std::bind(&fwdpy11::wf_rules::pick2, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, selfing_rate),
                std::bind(&fwdpy11::wf_rules::update, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, std::placeholders::_4,
                          std::placeholders::_5),
                mp);
            pop.N = N_next;
            um(pop.mutations, pop.fixations, pop.fixation_times,
               pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N);
            fitness.update(pop);
            auto wbar = rules.w(pop, fitness_callback);
			recorder(pop);
        }
}

// Evolve the population for some amount of time with mutation and
// recombination
void
evolve_singlepop_regions_cpp(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    py::array_t<std::uint32_t> popsizes, const double mu_neutral,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    fwdpy11::single_locus_fitness& fitness,
    fwdpy11::singlepop_temporal_sampler recorder, const double selfing_rate,
    const bool remove_selected_fixations = false)
{
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if (mu_neutral < 0.)
        {
            throw std::runtime_error("negative neutral mutation rate: "
                                     + std::to_string(mu_neutral));
        }
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    if (recrate < 0.)
        {
            throw std::runtime_error("negative recombination rate: "
                                     + std::to_string(recrate));
        }
    pop.mutations.reserve(std::ceil(
        std::log(2 * pop.N)
        * (4. * double(pop.N) * (mu_neutral + mu_selected)
           + 0.667 * (4. * double(pop.N) * (mu_neutral + mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    const auto mmodels = KTfwd::extensions::bind_dmm(
        mmodel, pop.mutations, pop.mut_lookup, rng.get(), mu_neutral,
        mu_selected, &pop.generation);
    ++pop.generation;
    auto rules = fwdpy11::wf_rules();
    if (remove_selected_fixations)
        {
            evolve_common(rng, pop, rules, popsizes, mu_neutral, mu_selected,
                          mmodels, recmap, fitness, recorder, selfing_rate,
                          std::true_type(),
                          fwdpy11::update_mutations_wrapper());
        }
    else
        {
            evolve_common(rng, pop, rules, popsizes, mu_neutral, mu_selected,
                          mmodels, recmap, fitness, recorder, selfing_rate,
                          KTfwd::remove_neutral(),
                          fwdpy11::update_mutations_n_wrapper());
        }
    --pop.generation;
}

PYBIND11_PLUGIN(wfevolve)
{
    py::module m("wfevolve", "example extending");

    m.def("evolve_singlepop_regions_cpp", &evolve_singlepop_regions_cpp);
    return m.ptr();
}
