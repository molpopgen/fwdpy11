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
}

PYBIND11_MODULE(wright_fisher_slocus_ts, m)
{
    m.doc() = "Evolution under a Wright-Fisher model using tree sequences.";

    m.def("WFSlocusPop_ts", &wfSlocusPop_ts);
}
