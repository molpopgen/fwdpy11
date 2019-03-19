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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/Mutation.hpp>

namespace py = pybind11;

struct flattened_Mutation
{
    double pos, s, h;
    fwdpp::uint_t g;
    decltype(fwdpy11::Mutation::xtra) label;
    std::int16_t neutral;
};

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);
PYBIND11_MAKE_OPAQUE(std::vector<flattened_Mutation>);

void
init_MutationVector(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(flattened_Mutation, pos, s, h, g, label, neutral);

    py::bind_vector<std::vector<fwdpy11::Mutation>>(
        m, "MutationVector",
        "C++ representation of a list of "
        ":class:`fwdpy11.Mutation`.  "
        "Typically, access will be read-only.",
        py::module_local(false))
        .def("array",
             [](const std::vector<fwdpy11::Mutation>& mc) {
                 std::vector<flattened_Mutation> rv;
                 rv.reserve(mc.size());
                 for (auto&& m : mc)
                     {
                         rv.push_back(flattened_Mutation{ m.pos, m.s, m.h, m.g,
                                                          m.xtra, m.neutral });
                     }
                 return rv;
             },
             R"delim(
        :rtype: :class:`fwdpy11.VecMutationDtype`.
        
        The return value should be coerced into a Numpy 
        array for processing.

        .. versionadded: 0.1.2
        )delim")
        .def(py::pickle(
            [](const std::vector<fwdpy11::Mutation>& mutations) {
                py::list rv;
                for (auto&& i : mutations)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpy11::Mutation> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpy11::Mutation>());
                    }
                return rv;
            }));

    py::bind_vector<std::vector<flattened_Mutation>>(
        m, "MutationDtypeVector", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the data fields in a "
        ":class:`fwdpy11.Mutation`.

        .. versionadded: 0.1.2
        )delim");
}
