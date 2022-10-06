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

std::vector<double>
convert_migmatrix(py::array_t<double> migmatrix)
{
    std::vector<double> migrates;
    auto r = migmatrix.unchecked<2>();
    for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
        {
            for (decltype(i) j = 0; j < r.shape(1); ++j)
                {
                    migrates.push_back(r(i, j));
                }
        }
    return migrates;
}

void
init_SetMigrationRate(py::module& m)
{
    py::class_<SMR>(m, "_ll_SetMigrationRates")
        .def(py::init<std::uint32_t, std::int32_t, std::vector<double>>(),
             py::arg("when"), py::arg("deme"), py::arg("migrates"))
        .def(py::init([](std::uint32_t when, py::array_t<double> migmatrix) {
                 auto m = convert_migmatrix(migmatrix);
                 return fwdpy11::discrete_demography::SetMigrationRates(when,
                                                                        std::move(m));
             }),
             py::arg("when"), py::arg("migmatrix"));
}
