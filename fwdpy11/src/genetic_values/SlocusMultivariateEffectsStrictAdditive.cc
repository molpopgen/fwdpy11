#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/SlocusMultivariateEffectsStrictAdditive.hpp>

namespace py = pybind11;

void
init_SlocusMultivariateEffectsStrictAdditive(py::module& m)
{
    py::class_<fwdpy11::SlocusMultivariateEffectsStrictAdditive,
               fwdpy11::SlocusPopMultivariateGeneticValueWithMapping>(
        m, "SlocusMultivariateEffectsStrictAdditive")
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::MultivariateGeneticValueToFitnessMap&>(),
             py::arg("ndimensions"), py::arg("focal_trait"),
             py::arg("genetic_values_to_fitness_map"))
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::MultivariateGeneticValueToFitnessMap&,
                      const fwdpy11::GeneticValueNoise&>(),
             py::arg("ndimensions"), py::arg("focal_trait"),
             py::arg("genetic_values_to_fitness_map"), py::arg("noise"))
        .def(py::pickle(
            [](const fwdpy11::SlocusMultivariateEffectsStrictAdditive& self) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    self.pickle(), p.attr("dumps")(self.gv2w->clone(), -1),
                    p.attr("dumps")(self.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid tuple size");
                    }
                auto data = t[0].cast<py::tuple>();
                auto ndim = data[0].cast<std::size_t>();
                auto focal_trait = data[1].cast<std::size_t>();
                auto p = py::module::import("pickle");
                auto gv2w = p.attr("loads")(t[1]);
                auto noise = p.attr("loads")(t[2]);
                return fwdpy11::SlocusMultivariateEffectsStrictAdditive(
                    ndim, focal_trait,
                    gv2w.cast<const fwdpy11::
                                  MultivariateGeneticValueToFitnessMap&>(),
                    noise.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
