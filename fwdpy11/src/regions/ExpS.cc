#include <pybind11/pybind11.h>
#include <fwdpy11/regions/ExpS.hpp>

namespace py = pybind11;

void
init_ExpS(py::module& m)
{
    py::class_<fwdpy11::ExpS, fwdpy11::Sregion>(m, "ExpS")
        .def(py::init<double, double, double, double, double, bool,
                      std::uint16_t, double>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("mean"), py::arg("h") = 1.0,
             py::arg("coupled") = true, py::arg("label") = 0,
             py::arg("scaling") = 1.0)
        .def_readonly("mean", &fwdpy11::ExpS::mean)
        .def_readonly("h", &fwdpy11::ExpS::dominance);
}

