#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/MultivariateGeneticValueToFitnessMap.hpp>

namespace py = pybind11;

void
init_MultivariateGeneticValueToFitnessMap(py::module &m)
{
    py::class_<fwdpy11::MultivariateGeneticValueToFitnessMap>(
        m, "MultivariateGeneticValueToFitnessMap",
        "ABC for classes representing mappings of multivariate traits to "
        "fitness.");
}
