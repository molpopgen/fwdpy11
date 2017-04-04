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
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <queue>
#include <cmath>
#include <stdexcept>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/rules/qtrait.hpp>

namespace py = pybind11;

//Generation time, optimum, VS, sigE
using env = std::tuple<KTfwd::uint_t, double, double, double>;

static const std::size_t GEN = 0;
static const std::size_t OPTIMUM=1;
static const std::size_t VS = 2;
static const std::size_t SIGE=3;

// Evolve the population for some amount of time with mutation and
// recombination
void
evolve_singlepop_regions_qtrait_cpp(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    py::array_t<std::uint32_t> popsizes, const double mu_neutral,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    fwdpy11::singlepop_fitness& fitness,
    fwdpy11::singlepop_temporal_sampler recorder, const double selfing_rate,
    const std::vector<env>& environmental_changes)
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
    if (environmental_changes.empty())
        {
            throw std::runtime_error("empty list of environmental changes.");
        }
    for (auto&& ei : environmental_changes)
        {
            if (std::get<VS>(ei) < 0.)
                {
                    throw std::runtime_error("VS must be >= 0.");
                }
            if (std::get<SIGE>(ei) < 0.)
                {
                    throw std::runtime_error(
                        "environmental noise must be >= 0.");
                }
        }
    const auto fitness_callback = fitness.callback();
    pop.mutations.reserve(std::ceil(
        std::log(2 * pop.N)
        * (4. * double(pop.N) * (mu_neutral + mu_selected)
           + 0.667 * (4. * double(pop.N) * (mu_neutral + mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    ++pop.generation;
    auto rules = fwdpy11::qtrait::qtrait_model_rules(0., 0., 0.);
    // generate FIFO queue of env changes
    std::queue<env> env_q(std::deque<env>(environmental_changes.begin(),
                                          environmental_changes.end()));
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            if (!env_q.empty())
                {
                    auto next_env = env_q.front();
                    if (pop.generation >= std::get<GEN>(next_env))
                        {
                            rules.optimum = std::get<OPTIMUM>(next_env);
                            rules.VS = std::get<VS>(next_env);
                            rules.sigE = std::get<SIGE>(next_env);
                            env_q.pop();
                        }
                }
            fitness.update(pop);
            const auto N_next = popsizes.at(generation);
            double wbar = KTfwd::experimental::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, N_next, mu_neutral + mu_selected,
                KTfwd::extensions::bind_dmm(
                    mmodel, pop.mutations, pop.mut_lookup, rng.get(),
                    mu_neutral, mu_selected, pop.generation),
                recmap, fitness_callback, pop.neutral, pop.selected,
                selfing_rate, rules);
            pop.N = N_next;
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
            recorder(pop);
        }
    --pop.generation;
}

PYBIND11_PLUGIN(wfevolve_qtrait)
{
    py::module m("wfevolve_qtrait", "example extending");

    m.def("evolve_singlepop_regions_qtrait_cpp",
          &evolve_singlepop_regions_qtrait_cpp);
    return m.ptr();
}
