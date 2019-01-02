#include <pybind11/pybind11.h>
#include <fwdpy11/regions/ConstantS.hpp>

namespace py = pybind11;

void
init_ConstantS(py::module& m)
{
    py::class_<fwdpy11::ConstantS, fwdpy11::Sregion>(m, "ConstantS")
        .def(py::init<double, double, double, double, double, bool,
                      std::uint16_t, double>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("s"), py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0)
        .def_readonly("esize", &fwdpy11::ConstantS::esize)
        .def_readonly("s", &fwdpy11::ConstantS::esize)
        .def_readonly("h", &fwdpy11::ConstantS::dominance);
}

