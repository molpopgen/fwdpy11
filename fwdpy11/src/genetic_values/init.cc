#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_MultivariateGeneticValueToFitnessMap(py::module&);
void init_MultivariateGSS(py::module&);
void init_MultivariateGSSmo(py::module&);
void init_DiploidPopulationMultivariateGeneticValueWithMapping(py::module&);
void init_DiploidMultivariateEffectsStrictAdditive(py::module&);

void
init_base_classes(py::module& m)
{
    init_MultivariateGeneticValueToFitnessMap(m);
    init_DiploidPopulationMultivariateGeneticValueWithMapping(m);
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
    init_DiploidMultivariateEffectsStrictAdditive(m);
}

void
init_genetic_values(py::module& m)
{
    init_base_classes(m);
    init_gvalue_to_fitness_classes(m);
    init_genetic_value_classes(m);
}

