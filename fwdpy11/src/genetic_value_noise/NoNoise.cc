#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/noise.hpp>

namespace py = pybind11;

void
init_NoNoise(py::module& m)
{
    py::class_<fwdpy11::NoNoise, fwdpy11::GeneticValueNoise>(
        m, "NoNoise", "Type implying no random effects on genetic values.")
        .def(py::init<>())
        .def(py::pickle(
            [](const fwdpy11::NoNoise& o) -> py::object { return o.pickle(); },
            [](py::object& o) { return fwdpy11::NoNoise::unpickle(o); }));
}
