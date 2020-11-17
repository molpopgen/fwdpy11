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
#include <cstdint>
#include <pybind11/pybind11.h>
#include <fwdpy11/discrete_demography/MassMigration.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

void
init_MassMigration(py::module& m)
{
    py::class_<ddemog::MassMigration>(m, "_ll_MassMigration")
        .def(py::init([](std::uint32_t when, std::int32_t source,
                         std::int32_t destination, double fraction,
                         bool move_individuals, bool resets_growth_rate) {
                 return ddemog::MassMigration(when, source, destination, 0, -1, fraction,
                                              move_individuals, false,
                                              resets_growth_rate);
             }),
             py::arg("when"), py::arg("source"), py::arg("destination"),
             py::arg("fraction"), py::arg("move_individuals"),
             py::arg("resets_growth_rate"));
}
