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
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fitness)
{
    py::module m("fitness", "Fitness models.");

    py::class_<fwdpy11::singlepop_fitness>(m, "SpopFitness",
                                           R"delim(
            A fitness function or trait value function
            for a single-deme, single-region simulation (
            :class:`fwdpy11.fwdpy11_types.Spop`).
            
            Current fitness ovjects derived from this type are:
            
            * :class:`fwdpy11.fitness.SpopAdditive`
            * :class:`fwdpy11.fitness.SpopMult`
            
            Current trait value objects derived from this type are:

            * :class:`fwdpy11.trait_values.SpopAdditiveTrait`
            * :class:`fwdpy11.trait_values.SpopMultTrait`
            * :class:`fwdpy11.trait_values.SpopGBRTrait`
            )delim");

    py::class_<fwdpy11::singlepop_mult_wrapper, fwdpy11::singlepop_fitness>(
        m, "SpopMult", R"delim(
        Multiplicative fitness for single-deme simulations.
        At a single mutation, fitness is 0, 1+sh, 1+scaling*s
        for genotypes AA, Aa, and aa, respectively. The scaling
        parameter is a constructor argument.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.SpopMult(1.0)
                       )delim")
        .def(py::init<double>(), py::arg("scaling"));

    py::class_<fwdpy11::singlepop_additive_wrapper,
               fwdpy11::singlepop_fitness>(m, "SpopAdditive", R"delim(
        Additive fitness for single-deme simulations.
        Fitness is max(0,1 + :math:`\sum_{i} x_i`), 
        where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
        for AA, Aa, and aa, respectively.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.SpopAdditive(2.0)
        )delim")
        .def(py::init<double>(), py::arg("scaling"));


    return m.ptr();
}
