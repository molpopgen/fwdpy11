#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_MultivariateGeneticValueToFitnessMap(py::module&);
void init_MultivariateGSS(py::module&);
void init_MultivariateGSSmo(py::module&);
void init_SlocusPopMultivariateGeneticValueWithMapping(py::module&);
void init_SlocusMultivariateEffectsStrictAdditive(py::module&);

void
init_genetic_values(py::module& m)
{
    init_MultivariateGeneticValueToFitnessMap(m);
    init_MultivariateGSS(m);
    init_MultivariateGSSmo(m);
    init_SlocusPopMultivariateGeneticValueWithMapping(m);
    init_SlocusMultivariateEffectsStrictAdditive(m);
}

