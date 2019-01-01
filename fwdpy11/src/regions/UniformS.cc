#include <pybind11/pybind11.h>
#include <fwdpy11/regions/UniformS.hpp>

namespace py = pybind11;

void
init_UniformS(py::module& m)
{
    py::class_<fwdpy11::UniformS, fwdpy11::Sregion>(m, "UniformS")
        .def(py::init<double, double, double, double, double, double, bool,
                      std::uint16_t, double>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("lo"),
             py::arg("hi"), py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0)
        .def_readonly("lo", &fwdpy11::UniformS::lo)
        .def_readonly("hi", &fwdpy11::UniformS::hi)
        .def_readonly("h", &fwdpy11::UniformS::dominance);
}

