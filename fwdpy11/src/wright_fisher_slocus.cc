// Wright-Fisher simulation for a fwdpy11::SlocusPop
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
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
// TODO: move concept of random noise to new set of headers.
#include <fwdpy11/evolve/qtrait_api.hpp>

namespace py = pybind11;

void
wfSlocusPop(
    const fwdpy11::GSLrng_t &rng, const fwdpy11::SlocusPop &pop,
    py::array_t<std::uint32_t> popsizes, const double mu_neutral,
    const double mu_selected, const double recrate,
    const fwdpp::extensions::discrete_mut_model<fwdpy11::SlocusPop::mcont_t>
        &mmodel,
    const fwdpp::extensions::discrete_rec_model &rmodel,
    fwdpy11::SlocusPopGeneticValue &genetic_value_fxn,
    fwdpy11::SlocusPopemporal_sampler recorder, const double selfing_rate,
    fwdpy11::single_locus_noise_function noise, py::object noise_updater,
    const bool remove_selected_fixations)
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
    const auto num_generations = popsizes.size();
    if (!num_generations)
        {
            throw std::invalid_argument("empty list of population sizes");
        }
}
