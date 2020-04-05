#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void
init_DiploidGeneticValue(py::module& m)
{
    py::class_<fwdpy11::DiploidGeneticValue>(
        m, "DiploidGeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.DiploidPopulation`")
        .def(
            "__call__",
            [](const fwdpy11::DiploidGeneticValue& gv, const std::size_t diploid_index,
               const fwdpy11::DiploidPopulation& pop) {
                return gv.calculate_gvalue(pop.diploid_metadata[diploid_index].label,
                                           pop.diploid_metadata[diploid_index], pop);
            },
            R"delim(
             :param diploid_index: The index of the individual to calculate.
             :type diploid_index: int >= 0
             :param pop: The population object containing the individual.
             :type pop: :class:`fwdpy11.DiploidPopulation`
             :return: The genetic value of an individual.
             :rtype: float
             )delim",
            py::arg("diploid_index"), py::arg("pop"))
        .def(
            "fitness",
            [](const fwdpy11::DiploidGeneticValue& gv, const std::size_t diploid_index,
               const fwdpy11::DiploidPopulation& pop) {
                return gv.genetic_value_to_fitness(pop.diploid_metadata[diploid_index]);
            },
            R"delim(
        :param diploid_index: The index of the individual
        :type diploid_index: int >= 0
        :param pop: The population containing the individual
        :type pop: :class:`fwdpy11.DiploidPopulation`
        :return: The fitness of an individual.
        :rtype: float
        )delim",
            py::arg("diploid_index"), py::arg("pop"))
        .def_property_readonly(
            "shape",
            [](const fwdpy11::DiploidGeneticValue& self) { return self.shape(); },
            R"delim(
                               Return the dimensions of the genetic values.

                               .. versionadded:: 0.3.0
                               )delim")
        .def_readonly("genetic_values", &fwdpy11::DiploidGeneticValue::gvalues,
                      R"delim(
                      Return the list of genetic values.

                      .. versionadded:: 0.3.0
                      )delim")
        .def_property_readonly(
            "gvalue_to_fitness",
            [](const fwdpy11::DiploidGeneticValue& o) { return o.gv2w->clone(); },
            "Access the genetic value to fitness map.")
        .def_property_readonly(
            "noise",
            [](const fwdpy11::DiploidGeneticValue& o) { return o.noise_fxn->clone(); },
            "Access the random noise funcion")
        .def_property_readonly(
            "maps_to_fitness",
            [](const fwdpy11::DiploidGeneticValue& self) {
                return self.gv2w->isfitness;
            },
            R"delim(
        Returns True if object represents a mapping directly to fitness, and
        False otherwise.

        .. versionadded:: 0.7.0
        )delim")
        .def_property_readonly(
            "maps_to_trait_value",
            [](const fwdpy11::DiploidGeneticValue& self) {
                return !self.gv2w->isfitness;
            },
            R"delim(
        Returns True if object represents a trait value, and
        False otherwise.

        .. versionadded:: 0.7.0
        )delim");
}
