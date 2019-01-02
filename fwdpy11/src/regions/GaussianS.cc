#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GaussianS.hpp>

namespace py = pybind11;

void
init_GaussianS(py::module& m)
{
    py::class_<fwdpy11::GaussianS, fwdpy11::Sregion>(m, "GaussianS")
        .def(py::init<double, double, double, double, double, bool,
                      std::uint16_t, double>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("sd"),
              py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0)
        .def_readonly("sd", &fwdpy11::GaussianS::sd)
        .def_readonly("h", &fwdpy11::GaussianS::dominance);
}

