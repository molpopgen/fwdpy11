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
#include <fwdpy11/genetic_value_to_fitness/Optimum.hpp>

namespace py = pybind11;
using namespace py::literals;

static const auto CLASS_DOCSTRING =
    R"delim(
Parameters for a trait optimum.

.. versionadded:: 0.7.1
)delim";

static const auto INIT_3 =
    R"delim(
:param when: The time when the optimum shifts
:type when: int
:param optimum: The trait value
:type optimum: float
:param VS: Strength of stabilizing selection
:type VS: float
)delim";

static const auto INIT_2 =
    R"delim(
:param optimum: The trait value
:type optimum: float
:param VS: Strength of stabilizing selection
:type VS: float
)delim";

void
init_Optimum(py::module& m)
{
    py::class_<fwdpy11::Optimum>(m, "Optimum", CLASS_DOCSTRING)
        .def(py::init<std::uint32_t, double, double>(), py::arg("when"),
             py::arg("optimum"), py::arg("VS"), INIT_3)
        .def(py::init<double, double>(), py::arg("optimum"), py::arg("VS"), INIT_2)
        .def_readonly("when", &fwdpy11::Optimum::when)
        .def_readonly("optimum", &fwdpy11::Optimum::opt)
        .def_readonly("VS", &fwdpy11::Optimum::VW)
        .def("__repr__",
             [](const fwdpy11::Optimum& self) {
                 return "Optimum(when={}, optimum={}, VS={})"_s.format(
                     self.when, self.opt, self.VW);
             })
        .def(py::pickle(
            [](const fwdpy11::Optimum& self) {
                return py::make_tuple(self.when, self.opt, self.VW);
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid tuple size");
                    }
                return fwdpy11::Optimum(t[0].cast<std::uint32_t>(), t[1].cast<double>(),
                                        t[2].cast<double>());
            }));
}
