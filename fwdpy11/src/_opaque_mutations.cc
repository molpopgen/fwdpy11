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
#include <fwdpp/sugar/popgenmut.hpp>

namespace py = pybind11;

struct flattened_popgenmut
{
    double pos, s, h;
    KTfwd::uint_t g;
    decltype(KTfwd::popgenmut::xtra) label;
    std::int8_t neutral;
};

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::popgenmut>);
PYBIND11_MAKE_OPAQUE(std::vector<flattened_popgenmut>);

PYBIND11_MODULE(_opaque_mutations, m)
{
    m.doc()
        = "Expose C++ vectors of Mutation objects to Python without copies.";

    PYBIND11_NUMPY_DTYPE(flattened_popgenmut, pos, s, h, g, label, neutral);

    py::bind_vector<std::vector<KTfwd::popgenmut>>(
        m, "VecMutation",
        "C++ representation of a list of "
        ":class:`fwdpy11.Mutation`.  "
        "Typically, access will be read-only.",
        py::module_local(false))
        .def("array",
             [](const std::vector<KTfwd::popgenmut>& mc) {
                 std::vector<flattened_popgenmut> rv;
                 rv.reserve(mc.size());
                 for (auto&& m : mc)
                     {
                         rv.push_back(flattened_popgenmut{
                             m.pos, m.s, m.h, m.g, m.xtra, m.neutral });
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
            [](const std::vector<KTfwd::popgenmut>& mutations) {
                py::list rv;
                for (auto&& i : mutations)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<KTfwd::popgenmut> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<KTfwd::popgenmut>());
                    }
                return rv;
            }));

    py::bind_vector<std::vector<flattened_popgenmut>>(
        m, "VecMutationDtype", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the data fields in a "
        ":class:`fwdpy11.Mutation`.

        .. versionadded: 0.1.2
        )delim");
}
