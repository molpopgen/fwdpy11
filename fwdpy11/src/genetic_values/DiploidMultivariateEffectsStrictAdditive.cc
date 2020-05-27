#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidMultivariateEffectsStrictAdditive.hpp>

namespace py = pybind11;

void
init_DiploidMultivariateEffectsStrictAdditive(py::module& m)
{
    py::class_<fwdpy11::DiploidMultivariateEffectsStrictAdditive,
               fwdpy11::DiploidGeneticValue>(m, "_ll_StrictAdditiveMultivariateEffects")
        .def(py::init([](std::size_t ndimensions, std::size_t focal_trait,
                         const fwdpy11::GeneticValueIsTrait& gvalue_to_fitness,
                         py::object noise) {
                 if (noise.is_none())
                     {
                         return fwdpy11::DiploidMultivariateEffectsStrictAdditive(
                             ndimensions, focal_trait, gvalue_to_fitness);
                     }
                 return fwdpy11::DiploidMultivariateEffectsStrictAdditive(
                     ndimensions, focal_trait, gvalue_to_fitness,
                     noise.cast<const fwdpy11::GeneticValueNoise&>());
             }),
             py::arg("ndimensions"), py::arg("focal_trait"),
             py::arg("gvalue_to_fitness"), py::arg("noise"));
}
