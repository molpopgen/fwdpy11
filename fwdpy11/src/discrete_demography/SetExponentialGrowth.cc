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
#include <sstream>
#include <pybind11/pybind11.h>
#include <fwdpy11/discrete_demography/SetExponentialGrowth.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;
using SEG = ddemog::SetExponentialGrowth;

static const auto INIT_DOCSTRING = R"delim(
:param when: The generation when the event occurs
:type when: int
:param deme: The deme whose growth rate will change
:type deme: int
:param G: The new growth rate
:type G: float
)delim";

void
init_SetExponentialGrowth(py::module& m)
{
    py::class_<SEG>(m, "SetExponentialGrowth",
                    R"delim(
        Set the growth rate of a deme at a given time.
        
        .. versionadded:: 0.5.3
        )delim")
        .def(py::init<decltype(SEG::when), decltype(SEG::deme),
                      decltype(SEG::G)>(),
             py::arg("when"), py::arg("deme"), py::arg("G"), INIT_DOCSTRING)
        .def_readonly("when", &SEG::when)
        .def_readonly("deme", &SEG::deme)
        .def_readonly("G", &SEG::G)
        .def("__repr__",
             [](const SEG& self) {
                 std::ostringstream o;
                 o << "SetExponentialGrowth(when=" << self.when
                   << ", deme=" << self.deme << ", G=" << self.G << ')';
                 return o.str();
             })
        .def(py::pickle(
            [](const SEG& self) {
                return py::make_tuple(self.when, self.deme, self.G);
            },
            [](py::tuple t) {
                return SEG(t[0].cast<decltype(SEG::when)>(),
                           t[1].cast<decltype(SEG::deme)>(),
                           t[2].cast<decltype(SEG::G)>());
            }));
}
