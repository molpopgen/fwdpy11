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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/forward_types.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::gamete>);

PYBIND11_MODULE(_opaque_gametes, m)
{
    m.doc() = "Expose C++ vectors of gametes to Python without copies.";

    py::bind_vector<std::vector<KTfwd::gamete>>(
        m, "VecGamete", py::module_local(false),
        "C++ representations of a list of "
        ":class:`fwdpy11.Gamete`.  "
        "Typically, access will be read-only.")
        .def(py::pickle(
            [](const std::vector<KTfwd::gamete>& gametes) {
                py::list rv;
                for (auto&& g : gametes)
                    {
                        rv.append(g);
                    }
                return rv;
            },
            [](const py::list l) {
                std::vector<KTfwd::gamete> rv;
                for (auto&& li : l)
                    {
                        rv.push_back(li.cast<KTfwd::gamete>());
                    }
                return rv;
            }));
}
