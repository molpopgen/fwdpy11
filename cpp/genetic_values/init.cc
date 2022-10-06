#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_DiploidGeneticValue(py::module&);
void init_PyDiploidGeneticValue(py::module&);
void init_Additive(py::module&);
void init_Multiplicative(py::module&);
void init_GBR(py::module&);
void init_DiploidMultivariateEffectsStrictAdditive(py::module&);
void init_dgvalue_pointer_vector(py::module&);

void
init_base_classes(py::module& m)
{
    init_DiploidGeneticValue(m);
    init_PyDiploidGeneticValue(m);
}

void
init_genetic_value_classes(py::module& m)
{
    init_Additive(m);
    init_Multiplicative(m);
    init_GBR(m);
    init_DiploidMultivariateEffectsStrictAdditive(m);
}

void
init_genetic_values(py::module& m)
{
    init_base_classes(m);
    init_genetic_value_classes(m);
    init_dgvalue_pointer_vector(m);
}

