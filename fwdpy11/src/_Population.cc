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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/Population.hpp>

namespace py = pybind11;

namespace
{
    static const auto MCOUNTS_DOCSTRING = R"delim(
    List of number of occurrences of elements in 
    a population objecst "mutations" container.

    The values are unsigned 32-bit integers.  

    .. note::
        Some values may be 0.  These represent *extinct* variants.  You will typically want to avoid processing such mutations.
)delim";

    static const auto FIXATIONS_DOCSTRING
        = R"delim(A :class:`fwdpy11.VecMutation` of fixed variants.)delim";

    static const auto FIXATION_TIMES_DOCSTRING =
        R"delim(A list of fixation times corresponding to the elements in "fixations" for this type.)delim";

    static const auto GAMETES_DOCSTRING
        = R"delim(A :class:`fwdpy11.VecGamete`.)delim";

    static const auto MUTATIONS_DOCSTRING = R"delim(
    List of :class:`fwdpy11.Mutation`.

    .. note:: 
        This list contains **both** extinct *and* extant mutations.  
        To distinguish them, use the locations of nonzero values in "mcounts" 
        for an instance of this type."
    )delim";

    static const auto POPDATA_DOCSTRING
        = "Python object that may be written to by a simulation. Any data "
          "written should be documented by the simulation function.\n\n.. "
          "versionadded:: 0.1.4";

    static const auto POPDATA_USER_DOCSTRING
        = "A Python object with read-write access.\n\n.. versionadded:: 0.1.4";
}

PYBIND11_MAKE_OPAQUE(fwdpy11::Population::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::mcont_t);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::uint_t>);

PYBIND11_MODULE(_Population, m)
{
    py::bind_vector<std::vector<fwdpp::uint_t>>(
        m, "VecUint32", "Vector of unsigned 32-bit integers.",
        py::buffer_protocol())
        .def(py::pickle(
            [](const std::vector<fwdpp::uint_t>& v) {
                py::list rv;
                for (auto&& i : v)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpp::uint_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpp::uint_t>());
                    }
                return rv;
            }));

    py::class_<fwdpy11::Population>(m, "_Population")
        .def_readonly("N", &fwdpy11::Population::N)
        .def_readonly("generation", &fwdpy11::Population::generation)
        .def_readonly("mutations", &fwdpy11::Population::mutations,
                      MUTATIONS_DOCSTRING)
        .def_readonly("mcounts", &fwdpy11::Population::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly("gametes", &fwdpy11::Population::gametes,
                      GAMETES_DOCSTRING)
        .def_readonly("fixations", &fwdpy11::Population::fixations,
                      FIXATIONS_DOCSTRING)
        .def_readonly("fixation_times", &fwdpy11::Population::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def_readonly("popdata", &fwdpy11::Population::popdata,
                      POPDATA_DOCSTRING)
        .def_readonly("popdata_user", &fwdpy11::Population::popdata_user,
                      POPDATA_USER_DOCSTRING);
}
