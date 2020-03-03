#include <fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueToFitnessMap(py::module& m)
{
    py::class_<fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueToFitnessMap",
        "ABC for functions translating genetic values into fitness.")
        .def_property_readonly("shape",&fwdpy11::GeneticValueToFitnessMap::shape);
}

