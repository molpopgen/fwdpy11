#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/poisson_point.hpp>

namespace py = pybind11;

void
init_PoissonPoint(py::module& m)
{
    py::class_<fwdpp::poisson_point, fwdpp::genetic_map_unit>(m, "_ll_PoissonPoint")
        .def(py::init<double, double>(), py::arg("position"), py::arg("mean"));
}

