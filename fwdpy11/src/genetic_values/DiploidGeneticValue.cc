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
        .def_property_readonly(
            "shape",
            [](const fwdpy11::DiploidGeneticValue& self) {
                if (self.total_dim > 1 && self.total_dim != self.gvalues.size())
                    {
                        throw std::runtime_error("dimensionality mismatch");
                    }
                return pybind11::make_tuple(self.total_dim);
            },
            R"delim(
                               Return the dimensions of the genetic values.

                               .. versionadded:: 0.3.0
                               )delim")
        .def_readonly("genetic_values", &fwdpy11::DiploidGeneticValue::gvalues,
                      R"delim(
                      Return the list of genetic values.

                      .. versionadded:: 0.3.0
                      )delim")
        //.def_property_readonly(
        //    "gvalue_to_fitness",
        //    [](const fwdpy11::DiploidGeneticValue& o) { return o.gv2w->clone(); },
        //    "Access the genetic value to fitness map.")
        //.def_property_readonly(
        //    "noise",
        //    [](const fwdpy11::DiploidGeneticValue& o) { return o.noise_fxn->clone(); },
        //    "Access the random noise funcion")
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
