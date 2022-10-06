#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_value_noise/NoNoise.hpp>

namespace py = pybind11;

void
init_NoNoise(py::module& m)
{
    py::class_<fwdpy11::NoNoise, fwdpy11::GeneticValueNoise>(m, "_ll_NoNoise")
        .def(py::init<>());
}
