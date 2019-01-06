#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/SlocusPopMultivariateGeneticValueWithMapping.hpp>

namespace py = pybind11;

void
init_SlocusPopMultivariateGeneticValueWithMapping(py::module& m)
{
    py::class_<fwdpy11::SlocusPopMultivariateGeneticValueWithMapping,
               fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusPopMultivariateGeneticValueWithMapping");
}
