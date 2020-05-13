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
#include <fwdpy11/discrete_demography/SetSelfingRate.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;
using SSR = ddemog::SetSelfingRate;

void
init_SetSelfingRate(py::module& m)
{
    py::class_<SSR>(m, "_ll_SetSelfingRate")
        .def(py::init<decltype(SSR::when), decltype(SSR::deme), decltype(SSR::S)>(),
             py::arg("when"), py::arg("deme"), py::arg("S"));
}
