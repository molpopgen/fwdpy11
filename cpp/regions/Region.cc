#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Region.hpp>

namespace py = pybind11;

void
init_Region(py::module& m)
{
    py::class_<fwdpy11::Region>(
        m, "_ll_Region")
        .def(py::init<double, double, double, bool, std::uint16_t>(), py::arg("beg"),
             py::arg("end"), py::arg("weight"), py::arg("coupled") = true,
             py::arg("label") = 0);
}

