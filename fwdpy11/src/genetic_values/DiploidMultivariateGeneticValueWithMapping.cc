#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidPopulationMultivariateGeneticValueWithMapping.hpp>

namespace py = pybind11;

void
init_DiploidPopulationMultivariateGeneticValueWithMapping(py::module& m)
{
    py::class_<fwdpy11::DiploidPopulationMultivariateGeneticValueWithMapping,
               fwdpy11::DiploidPopulationGeneticValue>(
        m, "MultivariateGeneticValueWithMapping",
        "ABC for multivariate traits.");
}
