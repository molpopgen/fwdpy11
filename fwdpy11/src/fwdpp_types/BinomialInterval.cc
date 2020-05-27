#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/binomial_interval.hpp>

namespace py = pybind11;

void
init_BinomialInterval(py::module& m)
{
    py::class_<fwdpp::binomial_interval, fwdpp::genetic_map_unit>(m,
                                                                  "_ll_BinomialInterval")
        .def(py::init<double, double, double>(), py::arg("beg"), py::arg("end"),
             py::arg("probability"));
}

