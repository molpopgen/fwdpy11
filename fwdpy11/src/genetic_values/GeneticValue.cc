#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValue(py::module& m)
{
    py::class_<fwdpy11::DiploidPopulationGeneticValue>(
        m, "GeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.DiploidPopulation`")
        .def("__call__",
             [](const fwdpy11::DiploidPopulationGeneticValue& gv,
                const std::size_t diploid_index,
                const fwdpy11::DiploidPopulation& pop) {
                 return gv.calculate_gvalue(diploid_index, pop);
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
        .def("fitness",
             [](const fwdpy11::DiploidPopulationGeneticValue& gv,
                const std::size_t diploid_index,
                const fwdpy11::DiploidPopulation& pop) {
                 return gv.genetic_value_to_fitness(
                     pop.diploid_metadata[diploid_index]);
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
            [](const fwdpy11::DiploidPopulationGeneticValue& self) {
                return self.shape();
            },
            R"delim(
                               Return the dimensions of the genetic values.

                               .. versionadded:: 0.3.0
                               )delim")
        .def_readonly("genetic_values",
                      &fwdpy11::DiploidPopulationGeneticValue::gvalues,
                      R"delim(
                      Return the list of genetic values.

                      .. versionadded:: 0.3.0
                      )delim");
}
