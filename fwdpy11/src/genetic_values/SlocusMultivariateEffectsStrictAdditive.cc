#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/SlocusMultivariateEffectsStrictAdditive.hpp>

namespace py = pybind11;

void
init_SlocusMultivariateEffectsStrictAdditive(py::module& m)
{
    py::class_<fwdpy11::SlocusMultivariateEffectsStrictAdditive,
               fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusMultivariateEffectsStrictAdditive")
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::GeneticValueIsTrait&>())
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::GeneticValueIsTrait&,
                      const fwdpy11::GeneticValueNoise&>());
}
