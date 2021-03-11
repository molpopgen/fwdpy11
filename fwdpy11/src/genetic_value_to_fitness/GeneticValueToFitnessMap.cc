#include <fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueToFitnessMap(py::module& m)
{
    py::class_<fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueToFitnessMap",
        "ABC for functions translating genetic values into fitness.")
        .def_property_readonly(
            "shape",
            [](const fwdpy11::GeneticValueToFitnessMap& self) {
                return pybind11::make_tuple(self.total_dim);
            },
            R"delim(
        Returns the shape (dimensonality) of the object

        .. versionadded:: 0.7.0
        )delim")
        .def_property_readonly(
            "maps_to_fitness",
            [](const fwdpy11::GeneticValueToFitnessMap& self) { return self.isfitness; },
            R"delim(
        Returns True if object represents a mapping directly to fitness, and
        False otherwise.

        .. versionadded:: 0.7.0
        )delim")
        .def_property_readonly(
            "maps_to_trait_value",
            [](const fwdpy11::GeneticValueToFitnessMap& self) {
                return !self.isfitness;
            },
            R"delim(
        Returns True if object represents a trait value, and
        False otherwise.

        .. versionadded:: 0.7.0
        )delim");
}

