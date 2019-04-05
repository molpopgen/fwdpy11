#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>

namespace py = pybind11;

void
init_Sregion(py::module& m)
{
    py::class_<fwdpy11::Sregion>(m, "Sregion",
                                 R"delim(
        Representation of a "region" in a simulation with a dominance term.

        This class is an ABC.

        .. versionchanged:: 0.13.a2
            Added "scaling" attribute.

        .. versionchanged:: 0.3.0
            Refactored from a pure Python class to a C++/pybind11 class
        )delim")
        .def_property_readonly(
            "b", [](const fwdpy11::Sregion& s) { return s.beg(); },
            "Beginning of region")
        .def_property_readonly(
            "e", [](const fwdpy11::Sregion& s) { return s.end(); },
            "End of region")
        .def_property_readonly(
            "w", [](const fwdpy11::Sregion& s) { return s.weight(); },
            "Weight")
        .def_property_readonly(
            "l", [](const fwdpy11::Sregion& s) { return s.label(); }, "Label")
        .def_property_readonly(
            "c", [](const fwdpy11::Sregion& s) { return s.region.coupled; },
            "Coupling parameter")
        .def_readonly("scaling", &fwdpy11::Sregion::scaling,
                      "Scaling parameter");
}

