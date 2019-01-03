#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_SlocusMultivariateEffectsStrictAdditive(py::module&);

void
init_genetic_values(py::module& m)
{
    init_SlocusMultivariateEffectsStrictAdditive(m);
}

