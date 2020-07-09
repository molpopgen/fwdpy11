#include <pybind11/pybind11.h>

namespace py = pybind11;

// Parameter classes
void init_Optimum(py::module&);
void init_PleiotropicOptima(py::module&);

// Univariate classes
void init_GeneticValueToFitnessMap(py::module&);
void init_GeneticValueIsTrait(py::module&);
void init_GeneticValueIsFitness(py::module&);
void init_GSSmo(py::module&);

// Multivariate classes
void init_MultivariateGSSmo(py::module&);

void
initialize_genetic_value_to_fitness(py::module& m)
{
    init_Optimum(m);
    init_PleiotropicOptima(m);

    init_GeneticValueToFitnessMap(m);
    init_GeneticValueIsTrait(m);
    init_GeneticValueIsFitness(m);
    init_GSSmo(m);

    init_MultivariateGSSmo(m);
}
