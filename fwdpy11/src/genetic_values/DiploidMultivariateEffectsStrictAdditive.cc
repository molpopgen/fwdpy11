#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidMultivariateEffectsStrictAdditive.hpp>

namespace py = pybind11;

void
init_DiploidMultivariateEffectsStrictAdditive(py::module& m)
{
    py::class_<fwdpy11::DiploidMultivariateEffectsStrictAdditive,
               fwdpy11::DiploidPopulationMultivariateGeneticValueWithMapping>(
        m, "StrictAdditiveMultivariateEffects",
        R"delim(
        Multivariate trait values under strictly additive effects.

        Calculate the trait value for a diploid in a :class:`fwdpy11.DiploidPopulation`
        for a multidimensional trait.

        This class is restricted to the case of simple additive effects, meaning
        that any dominance terms associated with mutations are ignored.

        During a simulation, :attr:`fwdpy11.DiploidMetadata.g` is filled with the 
        genetic value corresponding to a "focal" trait specified upon object construction.
        )delim")
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::MultivariateGeneticValueToFitnessMap&>(),
             py::arg("ndimensions"), py::arg("focal_trait"), py::arg("gv2w"),
             R"delim(
:param ndim: Number of trait dimensions
:type ndim: int
:param focal_trait: Index of the focal trait
:type focal_trait: int
:param gv2w: Function mapping trait value to fitness
:type gv2w: :class:`fwdpy11.MultivariateGeneticValueToFitnessMap`
            )delim")
        .def(py::init<std::size_t, std::size_t,
                      const fwdpy11::MultivariateGeneticValueToFitnessMap&,
                      const fwdpy11::GeneticValueNoise&>(),
             py::arg("ndimensions"), py::arg("focal_trait"),
             py::arg("genetic_values_to_fitness_map"), py::arg("noise"),
             R"delim(
:param ndim: Number of trait dimensions
:type ndim: int
:param focal_trait: Index of the focal trait
:type focal_trait: int
:param gv2w: Function mapping trait value to fitness
:type gv2w: :class:`fwdpy11.MultivariateGeneticValueToFitnessMap`
:param noise: Function adding random additive noise to trait value
:type noise: :class:`fwdpy11.GeneticValueNoise`
            )delim")

        .def(py::pickle(
            [](const fwdpy11::DiploidMultivariateEffectsStrictAdditive& self) {
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
                return fwdpy11::DiploidMultivariateEffectsStrictAdditive(
                    ndim, focal_trait,
                    gv2w.cast<const fwdpy11::
                                  MultivariateGeneticValueToFitnessMap&>(),
                    noise.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
