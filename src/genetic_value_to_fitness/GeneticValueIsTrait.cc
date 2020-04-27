#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueIsTrait(py::module &m)
{
    py::class_<fwdpy11::GeneticValueIsTrait,
               fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueIsTrait",
        "ABC for functions mapping genetic values representing traits to "
        "fitness.");
}
