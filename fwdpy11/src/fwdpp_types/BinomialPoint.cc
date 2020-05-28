#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/binomial_point.hpp>

namespace py = pybind11;
using namespace py::literals;

void
init_BinomialPoint(py::module& m)
{
    py::class_<fwdpp::binomial_point, fwdpp::genetic_map_unit>(m, "_ll_BinomialPoint")
        .def(py::init<double, double>(), py::arg("position"), py::arg("probability"));
}

