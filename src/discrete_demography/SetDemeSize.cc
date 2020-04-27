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
#include <fwdpy11/discrete_demography/SetDemeSize.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;
using SDS = ddemog::SetDemeSize;

static const auto INIT_DOCSTRING = R"delim(
:param when: The generation when the event occurs
:type when: int
:param deme: The deme whose size will change
:type deme: int
:param new_size: The new size
:type new_size: int
:param resets_growth_rate: (True) If deme size change resets growth rate to :data:`fwdpy11.NOGROWTH`
:type resets_growth_rate: bool
)delim";

void
init_SetDemeSize(py::module &m)
{
    py::class_<SDS>(m, "SetDemeSize",
                    R"delim(
        Set the size of a deme at a given time.
        
        .. versionadded:: 0.5.3
        )delim")
        .def(py::init<decltype(SDS::when), decltype(SDS::deme),
                      decltype(SDS::new_size),
                      decltype(SDS::resets_growth_rate)>(),
             py::arg("when"), py::arg("deme"), py::arg("new_size"),
             py::arg("resets_growth_rate") = true, INIT_DOCSTRING)
        .def_readonly("when", &SDS::when)
        .def_readonly("deme", &SDS::deme)
        .def_readonly("new_size", &SDS::new_size)
        .def_readonly("resets_growth_rate", &SDS::resets_growth_rate)
        .def("__repr__",
             [](const SDS &self) {
                 std::ostringstream o;
                 o << "SetDemeSize(when=" << self.when
                   << ", deme=" << self.deme << ", new_size=" << self.new_size
                   << ", resets_growth_rate=" << self.resets_growth_rate
                   << ")";
                 return o.str();
             })
        .def(py::pickle(
            [](const SDS &self) {
                return py::make_tuple(self.when, self.deme, self.new_size,
                                      self.resets_growth_rate);
            },
            [](py::tuple t) {
                return SDS(t[0].cast<decltype(SDS::when)>(),
                           t[1].cast<decltype(SDS::deme)>(),
                           t[2].cast<decltype(SDS::new_size)>(),
                           t[3].cast<decltype(SDS::resets_growth_rate)>());
            }));
}
