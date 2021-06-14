//
// Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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
#include "add_mutation.hpp"

namespace py = pybind11;

void
init_add_mutation(py::module& m)
{
    m.def("_add_mutation", &add_mutation);

    py::class_<new_mutation_data>(m, "ll_NewMutationData")
        .def(py::init<double, double, std::vector<double>, std::vector<double>,
                      decltype(fwdpy11::Mutation::xtra)>());
}

