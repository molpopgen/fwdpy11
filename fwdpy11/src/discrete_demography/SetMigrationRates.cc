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

#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpy11/discrete_demography/SetMigrationRates.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

using SMR = ddemog::SetMigrationRates;

static const auto INIT_DOCSTRING = R"delim(
:param when: The generation when the event occurs
:type when: int
:param deme: The row index of the migration matrix
:type when: int
:param migrates: The migration rates from `deme` to all other populations.
:type migrates: list
)delim";

static const auto INIT_DOCSTRING_NUMPY = R"delim(
:param when: The generation when the event occurs
:type when: int
:param migmatrix: A square matrix representing the migration matrix
:type migmatrix: numpy.ndarray
)delim";

void
init_SetMigrationRate(py::module& m)
{
    py::class_<SMR>(m, "SetMigrationRates",
                    R"delim(
        Set the migration parameters of a simulation at a given time.
        May be used to set either the migration rates from a given
        deme or the entire migration matrix.
        
        .. versionadded:: 0.5.3
        )delim")

        .def(py::init<std::uint32_t, std::int32_t, std::vector<double>>(),
             py::arg("when"), py::arg("deme"), py::arg("migrates"),
             INIT_DOCSTRING)
        .def(py::init<std::uint32_t, py::array_t<double>>(), py::arg("when"),
             py::arg("migmatrix"), INIT_DOCSTRING_NUMPY)
        .def_readonly("when", &SMR::when)
        .def_readonly("deme", &SMR::deme)
        .def_property_readonly(
            "migrates",
            [](const SMR& self) {
                if (self.deme == ddemog::NULLDEME)
                    {
                        auto dim = std::sqrt(self.migrates.size());
                        return fwdpy11::make_2d_ndarray_readonly(self.migrates,
                                                                 dim, dim);
                    }
                else
                    {
                        return fwdpy11::make_1d_ndarray_readonly(
                            self.migrates);
                    }
            })
        .def(py::pickle(
            [](const SMR& self) {
                return py::make_tuple(self.when, self.deme, self.migrates);
            },
            [](py::tuple t) {
                return SMR(t[0].cast<decltype(SMR::when)>(),
                           t[1].cast<decltype(SMR::deme)>(),
                           t[2].cast<decltype(SMR::migrates)>());
            }));
}
