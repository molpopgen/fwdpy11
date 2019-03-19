#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_GeneticValueToFitnessMap(py::module&);
void init_GeneticValueIsTrait(py::module&);
void init_GeneticValueIsFitness(py::module&);
void init_GSS(py::module&);
void init_GSSmo(py::module&);

// Multivariate classes
void init_MultivariateGeneticValueToFitnessMap(py::module&);
void init_MultivariateGSS(py::module&);
void init_MultivariateGSSmo(py::module&);

void
initialize_genetic_value_to_fitness(py::module& m)
{
    init_GeneticValueToFitnessMap(m);
    init_GeneticValueIsTrait(m);
    init_GeneticValueIsFitness(m);
    init_GSS(m);
    init_GSSmo(m);

    init_MultivariateGeneticValueToFitnessMap(m);
    init_MultivariateGSS(m);
    init_MultivariateGSSmo(m);
}
