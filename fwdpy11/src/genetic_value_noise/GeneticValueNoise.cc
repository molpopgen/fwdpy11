#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/noise.hpp>

namespace py = pybind11;

void
init_GeneticValueNoise(py::module& m)
{
    py::class_<fwdpy11::GeneticValueNoise>(
        m, "GeneticValueNoise",
        "ABC for noise classes affecting :class:`fwdpy11.DiploidPopulation`.");
}
