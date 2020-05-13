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
#include <fwdpy11/discrete_demography/MigrationMatrix.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

void
init_MigrationMatrix(py::module &m)
{
    py::class_<ddemog::MigrationMatrix>(m, "_ll_MigrationMatrix")
        .def(py::init([](py::array_t<double> m, const bool scaled) {
                 auto r = m.unchecked<2>();
                 if (r.shape(0) != r.shape(1))
                     {
                         throw std::invalid_argument("MigrationMatrix must be square");
                     }
                 auto d = r.data(0, 0);
                 std::vector<double> M(d, d + r.shape(0) * r.shape(1));
                 return ddemog::MigrationMatrix(std::move(M), r.shape(0), scaled);
             }),
             py::arg("migmatrix"), py::arg("scaled") = false);
}

