#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_GeneticValueNoise(py::module&);
void init_NoNoise(py::module&);
void init_GaussianNoise(py::module&);

void
initialize_genetic_value_noise(py::module& m)
{
    init_GeneticValueNoise(m);
    init_NoNoise(m);
    init_GaussianNoise(m);
}
