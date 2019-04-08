#include <fwdpy11/genetic_values/DiploidPopulationGeneticValueWithMapping.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueWithMapping(py::module& m)
{
    py::class_<fwdpy11::DiploidPopulationGeneticValueWithMapping,
               fwdpy11::DiploidPopulationGeneticValue>(
        m, "GeneticValueWithMapping",
        "ABC for genetic value calculations with flexible mapping of genetic value to "
        "fitness.")
        .def_property_readonly(
            "gvalue_to_fitness",
            [](const fwdpy11::DiploidPopulationGeneticValueWithMapping& o) {
                return o.gv2w->clone();
            },
            "Access the genetic value to fitness map.")
        .def_property_readonly(
            "noise",
            [](const fwdpy11::DiploidPopulationGeneticValueWithMapping& o) {
                return o.noise_fxn->clone();
            },
            "Access the random noise funcion");
}
