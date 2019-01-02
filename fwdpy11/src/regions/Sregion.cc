#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>

namespace py = pybind11;

void
init_Sregion(py::module& m)
{
    py::class_<fwdpy11::Sregion>(m, "Sregion")
        .def_property_readonly(
            "beg", [](const fwdpy11::Sregion& s) { return s.beg(); })
        .def_property_readonly(
            "end", [](const fwdpy11::Sregion& s) { return s.end(); })
        .def_property_readonly(
            "weight", [](const fwdpy11::Sregion& s) { return s.weight(); })
        .def_property_readonly(
            "label", [](const fwdpy11::Sregion& s) { return s.label(); })
        .def_readonly("scaling", &fwdpy11::Sregion::scaling);
}

