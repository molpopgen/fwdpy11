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

void
init_Optimum(py::module& m)
{
    py::class_<fwdpy11::Optimum>(m, "_ll_Optimum")
        .def(py::init([](double optimum, double VS, py::object when) {
                 if (when.is_none())
                     {
                         return fwdpy11::Optimum(optimum, VS);
                     }
                 return fwdpy11::Optimum(when.cast<std::uint32_t>(), optimum, VS);
             }),
             py::arg("optimum"), py::arg("VS"), py::arg("when"));
}
