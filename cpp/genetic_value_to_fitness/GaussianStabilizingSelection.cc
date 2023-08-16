#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/GaussianStabilizingSelection.hpp>

namespace py = pybind11;

void
init_GaussianStabilizingSelection(py::module& m)
{
    py::class_<fwdpy11::GaussianStabilizingSelection, fwdpy11::GeneticValueIsTrait>(
        m, "_ll_GaussianStabilizingSelection")
        .def(py::init([](const fwdpy11::GSSmo& single_trait) {
            return fwdpy11::GaussianStabilizingSelection(single_trait);
        }))
        .def(py::init([](const fwdpy11::MultivariateGSSmo& pleiotropy) {
            return fwdpy11::GaussianStabilizingSelection(pleiotropy);
        }));
}
