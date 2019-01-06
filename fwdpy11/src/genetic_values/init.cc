#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_MultivariateGeneticValueToFitnessMap(py::module&);
void init_MultivariateGSS(py::module&);
void init_MultivariateGSSmo(py::module&);
void init_SlocusPopMultivariateGeneticValueWithMapping(py::module&);
void init_SlocusMultivariateEffectsStrictAdditive(py::module&);

void
init_base_classes(py::module& m)
{
    init_MultivariateGeneticValueToFitnessMap(m);
    init_SlocusPopMultivariateGeneticValueWithMapping(m);
}

void
init_gvalue_to_fitness_classes(py::module& m)
{
    init_MultivariateGSS(m);
    init_MultivariateGSSmo(m);
}

void
init_genetic_value_classes(py::module& m)
{
    init_SlocusMultivariateEffectsStrictAdditive(m);
}

void
init_genetic_values(py::module& m)
{
    init_base_classes(m);
    init_gvalue_to_fitness_classes(m);
    init_genetic_value_classes(m);
}

