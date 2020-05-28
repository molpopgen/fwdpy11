#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/poisson_interval.hpp>

namespace py = pybind11;

void
init_PoissonInterval(py::module& m)
{
    py::class_<fwdpp::poisson_interval, fwdpp::genetic_map_unit>(m, "_ll_PoissonInterval")
        .def(py::init<double, double, double>(), py::arg("beg"), py::arg("end"),
             py::arg("mean"));
}

