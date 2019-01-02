#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>

namespace py = pybind11;

void
init_Sregion(py::module& m)
{
    py::class_<fwdpy11::Sregion>(m, "Sregion")
        .def_property_readonly(
            "b", [](const fwdpy11::Sregion& s) { return s.beg(); })
        .def_property_readonly(
            "e", [](const fwdpy11::Sregion& s) { return s.end(); })
        .def_property_readonly(
            "w", [](const fwdpy11::Sregion& s) { return s.weight(); })
        .def_property_readonly(
            "l", [](const fwdpy11::Sregion& s) { return s.label(); })
        .def_property_readonly(
            "c", [](const fwdpy11::Sregion& s) { return s.region.coupled; })
        .def_readonly("scaling", &fwdpy11::Sregion::scaling);
}

