#include <fwdpy11/genetic_value_to_fitness/GSS.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GSS(py::module& m)
{
    py::class_<fwdpy11::GSS, fwdpy11::GeneticValueIsTrait>(m, "_ll_GSS")
        .def(py::init<double, double>(), py::arg("optimum"), py::arg("VS"))
        .def(py::init<fwdpy11::Optimum>(), py::arg("optimum"));
}
