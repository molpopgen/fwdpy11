#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/SlocusPopMultivariateEffectsStrictAdditive.hpp>

namespace py = pybind11;

void
init_SlocusPopMultivariateEffectsStrictAdditive(py::module& m)
{
    py::class_<fwdpy11::SlocusPopMultivariateEffectsStrictAdditive,
               fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusPopMultivariateEffectsStrictAdditive");
}
