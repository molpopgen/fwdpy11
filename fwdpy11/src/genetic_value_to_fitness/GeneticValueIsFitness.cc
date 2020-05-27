#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsFitness.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueIsFitness(py::module& m)
{
    py::class_<fwdpy11::GeneticValueIsFitness, fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueIsFitness", "Type implying the the genetic value is fitness.")
        .def(py::init<std::size_t>())
        .def(py::pickle(
            [](const fwdpy11::GeneticValueIsFitness& g) {
                return py::make_tuple(g.total_dim);
            },
            [](py::tuple t) {
                std::size_t ndim = t[0].cast<std::size_t>();
                return fwdpy11::GeneticValueIsFitness(ndim);
            }));
}
