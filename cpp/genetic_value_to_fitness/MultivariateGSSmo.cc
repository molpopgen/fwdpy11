#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSSmo.hpp>

namespace py = pybind11;

void
init_MultivariateGSSmo(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "_ll_MultivariateGSSmo")
        .def(py::init<const std::vector<fwdpy11::PleiotropicOptima>&>(),
             py::arg("optima"));
}
