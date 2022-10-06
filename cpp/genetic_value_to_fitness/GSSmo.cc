#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/GSSmo.hpp>

namespace py = pybind11;

void
init_GSSmo(py::module& m)
{
    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueIsTrait>(m, "_ll_GSSmo")
        .def(py::init<std::vector<fwdpy11::Optimum>>(), py::arg("optima"));
}
