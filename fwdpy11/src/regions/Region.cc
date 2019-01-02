#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Region.hpp>

namespace py = pybind11;

void
init_Region(py::module& m)
{
    py::class_<fwdpy11::Region>(m, "Region")
        .def(py::init<double, double, double, bool, std::uint16_t>(), py::arg("beg"),
             py::arg("end"), py::arg("weight"), py::arg("coupled") = true,py::arg("label")=0)
        .def_readonly("beg", &fwdpy11::Region::beg)
        .def_readonly("end", &fwdpy11::Region::end)
        .def_readonly("weight", &fwdpy11::Region::weight)
        .def_readonly("coupled", &fwdpy11::Region::coupled)
        .def_readonly("label", &fwdpy11::Region::label);
}

