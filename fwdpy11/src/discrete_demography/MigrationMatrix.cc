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

static const auto INIT_DOCSTRING = R"delim(

:param migmatrix: A square matrix of non-negative floats.
:type migmatrix: numpy.ndarray
:param scaled: (True) If entries in `migmatrix` will be multiplied by deme sizes during simulation
:type scaled: bool
)delim";

void
init_MigrationMatrix(py::module &m)
{
    py::class_<ddemog::MigrationMatrix>(m, "MigrationMatrix",
                                         R"delim(
        The forward migration matrix for a simulation.

        .. versionadded:: 0.5.3
        )delim")
        .def(py::init([](py::array_t<double> m, const bool scaled) {
                 auto r = m.unchecked<2>();
                 if (r.shape(0) != r.shape(1))
                     {
                         throw std::invalid_argument(
                             "MigrationMatrix must be square");
                     }
                 auto d = r.data(0, 0);
                 std::vector<double> M(d, d + r.shape(0) * r.shape(1));
                 return ddemog::MigrationMatrix(std::move(M), r.shape(0),
                                                 scaled);
             }),
             py::arg("migmatrix"), py::arg("scale_during_simulation") = true,
             INIT_DOCSTRING)
        .def_property_readonly("shape",
                               [](ddemog::MigrationMatrix &self) {
                                   return pybind11::make_tuple(self.npops,
                                                               self.npops);
                               })
        .def_property_readonly("M",
                               [](const ddemog::MigrationMatrix &self) {
                                   auto M(self.M);
                                   return fwdpy11::make_2d_array_with_capsule(
                                       std::move(M), self.npops, self.npops);
                               })
        .def_readonly("scaled", &ddemog::MigrationMatrix::scaled)
        .def("_set_migration_rates",
             &ddemog::MigrationMatrix::set_migration_rates, py::arg("source"),
             py::arg("rates"))
        .def("_set_migration_rates",
             [](ddemog::MigrationMatrix &self, py::array_t<double> migrates) {
                 auto r = migrates.unchecked<2>();
                 std::vector<double> m;
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         for (decltype(i) j = 0; j < r.shape(1); ++j)
                             {
                                 m.push_back(r(i, j));
                             }
                     }
                 self.set_migration_rates(ddemog::NULLDEME, m);
             })
        .def(py::pickle(
            [](const ddemog::MigrationMatrix &self) {
                return py::make_tuple(self.M, self.npops, self.scaled);
            },
            [](pybind11::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid tuple size");
                    }
                auto M = t[0].cast<std::vector<double>>();
                auto n = t[1].cast<std::size_t>();
                auto s = t[2].cast<bool>();
                return ddemog::MigrationMatrix(std::move(M), n, s);
            }));
}

