#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_GeneticValue(py::module&);
void init_GeneticValueWithMapping(py::module&);
void init_Additive(py::module&);
void init_Multiplicative(py::module&);
void init_GBR(py::module&);
void init_DiploidPopulationMultivariateGeneticValueWithMapping(py::module&);
void init_DiploidMultivariateEffectsStrictAdditive(py::module&);
void init_dgvalue_pointer_vector(py::module&);

void
init_base_classes(py::module& m)
{
    init_GeneticValue(m);
    init_GeneticValueWithMapping(m);
    init_DiploidPopulationMultivariateGeneticValueWithMapping(m);
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

