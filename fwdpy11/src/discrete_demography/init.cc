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
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <fwdpy11/discrete_demography/constants.hpp>

namespace py = pybind11;

void init_discrete_demography_exceptions(py::module&);
void init_MigrationMatrix(py::module&);

void init_MassMigration(py::module&);
void init_SetDemeSize(py::module&);
void init_SetExponentialGrowth(py::module&);
void init_SetSelfingRate(py::module&);
void init_SetMigrationRate(py::module&);
void init_DiscreteDemography(py::module&);

void
init_discrete_demography(py::module& m)
{
    m.attr("NOGROWTH") = py::float_(fwdpy11::discrete_demography::NOGROWTH);

    // Concrete classes & functions now
    init_discrete_demography_exceptions(m);
    init_MigrationMatrix(m);

    init_MassMigration(m);
    init_SetDemeSize(m);
    init_SetExponentialGrowth(m);
    init_SetSelfingRate(m);
    init_SetMigrationRate(m);
    init_DiscreteDemography(m);
}
