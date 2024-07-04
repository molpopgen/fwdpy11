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

#include <functional>
#include <stdexcept>

#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/genetic_values/DiploidAdditive.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_Additive(py::module& m)
{
    py::class_<fwdpy11::DiploidAdditive, fwdpy11::DiploidGeneticValue>(m, "_ll_Additive")
        .def(py::init([](double scaling,
                         const fwdpy11::GeneticValueIsTrait* gvalue_to_fitness,
                         const fwdpy11::GeneticValueNoise* noise, std::size_t ndemes) {
                 if (gvalue_to_fitness != nullptr)
                     {
                         return fwdpy11::additive_trait_model(ndemes, scaling,
                                                              gvalue_to_fitness, noise);
                     }
                 return fwdpy11::additive_fitness_model(ndemes, scaling, noise);
             }),
             py::arg("scaling"), py::arg("gvalue_to_fitness"), py::arg("noise"),
             py::arg("ndemes"));
}
