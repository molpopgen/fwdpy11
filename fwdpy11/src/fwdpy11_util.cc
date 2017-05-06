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
#include "fwdpy11_util_add_mutation.hpp"

namespace py = pybind11;

PYBIND11_PLUGIN(util)
{
    py::module m("util", "Miscellaneous utilities for simulations.");

    m.def("add_mutation",
          [](const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
             const KTfwd::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label) {
              return add_mutation(rng, pop, ncopies, pos_s_h, label);
          },
          py::arg("rng"), py::arg("pop"), py::arg("ncopies"),
          py::arg("pos_esize_h"), py::arg("label") = 0,
          R"delim(
          Add a new mutation to a population.

          :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`.
          :param pop: A :class:`fwdpy11.fwdpy11_types.SlocusPop`.
          :param ncopies: The number of copies of the mutation.
          :param pos_esize_h: A tuple containing (pos,s,h) for the mutation.
          :param label: (0) An integer label for the mutation.

          .. note::
            The genotype frequencies for the new variant will
            be multinomially distributed around expected 
            Hardy-Weinberg proportions.

          .. note:: 
            An exception will occur if the new mutation 
            position is at a pre-existing position in 
            the population.
          )delim");

    m.def("add_mutation",
          [](const fwdpy11::GSLrng_t& rng, fwdpy11::multilocus_t& pop,
             const std::size_t locus, const KTfwd::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label) {
              return add_mutation(rng, pop, locus, ncopies, pos_s_h, label);
          },
          py::arg("rng"), py::arg("pop"), py::arg("locus"), py::arg("ncopies"),
          py::arg("pos_esize_h"), py::arg("label") = 0,
          R"delim(
          Add a new mutation to a population.

          :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`.
          :param pop: A :class:`fwdpy11.fwdpy11_types.MlocusPop`.
          :param locus: The locus in which the mutation will be located.
          :param ncopies: The number of copies of the mutation.
          :param pos_esize_h: A tuple containing (pos,s,h) for the mutation.
          :param label: (0) An integer label for the mutation.

          .. note::
            The genotype frequencies for the new variant will
            be multinomially distributed around expected 
            Hardy-Weinberg proportions.

          .. note:: 
            An exception will occur if the new mutation 
            position is at a pre-existing position in 
            the population.

          .. todo::
            Enforce that the new mutation position is 
            correct given pop.locus_boundaries
          )delim");


    return m.ptr();
}
