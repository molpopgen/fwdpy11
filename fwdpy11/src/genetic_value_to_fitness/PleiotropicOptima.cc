//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/genetic_value_to_fitness/PleiotropicOptima.hpp>

namespace py = pybind11;
using namespace py::literals;

static const auto INIT_3 =
    R"delim(
:param when: When the optima shift
:type when: int
:param optima: the new optima values
:type optima: list
:param VS: strength of stablizing selection
:type VS: float
)delim";

const auto INIT_2 =
    R"delim(
:param optima: the new optima values
:type optima: list
:param VS: strength of stablizing selection
:type VS: float
)delim";

void
init_PleiotropicOptima(py::module& m)
{
    py::class_<fwdpy11::PleiotropicOptima>(m, "PleiotropicOptima")
        .def(py::init<std::uint32_t, std::vector<double>, double>(), py::arg("when"),
             py::arg("optima"), py::arg("VS"), INIT_3)
        .def(py::init<std::vector<double>, double>(), py::arg("optima"), py::arg("VS"),
             INIT_2)
        .def_readonly("when", &fwdpy11::PleiotropicOptima::when)
        .def_readonly("optima", &fwdpy11::PleiotropicOptima::optima)
        .def_readonly("VS", &fwdpy11::PleiotropicOptima::VW)
        .def("__repr__",
             [](const fwdpy11::PleiotropicOptima& self) {
                 return "PleiotropicOptima(when={}, optima={}, VS={})"_s.format(
                     self.when, self.optima, self.VW);
             })
        .def(py::pickle(
            [](const fwdpy11::PleiotropicOptima& self) {
                return py::make_tuple(self.when, self.optima, self.VW);
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid tuple size");
                    }
                return fwdpy11::PleiotropicOptima(t[0].cast<double>(),
                                                  t[1].cast<std::vector<double>>(),
                                                  t[2].cast<double>());
            }));
}
