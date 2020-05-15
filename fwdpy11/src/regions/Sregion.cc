#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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

        .. versionchanged:: 0.8.0

            C++ version of class API minimized in favor of attrs.
        )delim")
        .def_property_readonly("dominance", &fwdpy11::Sregion::get_dominance,
                               "Dominance values.  Added in 0.7.0")
        .def_property_readonly("shape", &fwdpy11::Sregion::shape,
                               "Return shape.  Added in 0.7.0");
}

