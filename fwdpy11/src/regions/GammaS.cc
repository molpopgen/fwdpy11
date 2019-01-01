#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GammaS.hpp>

namespace py = pybind11;

void
init_GammaS(py::module& m)
{
    py::class_<fwdpy11::GammaS, fwdpy11::Sregion>(m, "GammaS")
        .def(py::init<double, double, double, double, double, double, bool,
                      std::uint16_t, double>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("mean"), py::arg("shape"), py::arg("h") = 1.0,
             py::arg("coupled") = true, py::arg("label") = 0,
             py::arg("scaling") = 1.0)
        .def_readonly("mean", &fwdpy11::GammaS::mean)
        .def_readonly("shape", &fwdpy11::GammaS::shape)
        .def_readonly("h", &fwdpy11::GammaS::dominance);
}
