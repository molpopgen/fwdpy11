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
#include <fwdpy11/types/Diploid.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidGenotype>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<fwdpy11::DiploidGenotype>>);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidMetadata>);

PYBIND11_MODULE(_opaque_diploids, m)
{
    m.doc() = "Expose C++ containers of diploids to Python without copies.";
    PYBIND11_NUMPY_DTYPE(fwdpy11::DiploidGenotype, first, second);

    py::bind_vector<fwdpy11::dipvector_t>(
        m, "VecDiploid", py::buffer_protocol(), py::module_local(false),
        "C++ representation of a list of "
        ":class:`fwdpy11."
        "DiploidGenotype`.  Typically, access will be read-only.")
        .def(py::pickle(
            [](const std::vector<fwdpy11::DiploidGenotype>& v) -> py::list {
                py::list rv;
                for (auto&& vi : v)
                    {
                        rv.append(vi);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpy11::DiploidGenotype> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpy11::DiploidGenotype>());
                    }
                return rv;
            }));

    py::bind_vector<std::vector<fwdpy11::dipvector_t>>(
        m, "VecVecDiploid", py::module_local(false),
        "Vector of "
        ":class:`fwdpy11.DiploidGenotype`.")
        .def(py::pickle(
            [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                py::list rv;
                for (auto&& i : diploids)
                    {
                        rv.append(i);
                    };
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpy11::dipvector_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpy11::dipvector_t>());
                    }
                return rv;
            }));

    PYBIND11_NUMPY_DTYPE(fwdpy11::DiploidMetadata, g, e, w, geography, label,
                         parents, deme, sex);

    py::bind_vector<std::vector<fwdpy11::DiploidMetadata>>(
        m, "VecDiploidMetaData", py::module_local(false),
        py::buffer_protocol(),
        R"delim(
        Container of diploid metadata.
        )delim");
}
