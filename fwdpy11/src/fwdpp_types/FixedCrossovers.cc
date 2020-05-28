#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/fixed_number_crossovers.hpp>

namespace py = pybind11;

void
init_FixedCrossovers(py::module& m)
{
    py::class_<fwdpp::fixed_number_crossovers, fwdpp::genetic_map_unit>(
        m, "_ll_FixedCrossovers")
        .def(py::init<double, double, int>(), py::arg("beg"), py::arg("end"),
             py::arg("num_xovers"));
}

