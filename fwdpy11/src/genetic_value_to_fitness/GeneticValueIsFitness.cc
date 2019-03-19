#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_GeneticValueIsFitness(py::module& m)
{
    py::class_<fwdpy11::GeneticValueIsFitness,
               fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueIsFitness",
        "Type implying the the genetic value is fitness.")
        .def(py::init<>())
        .def(py::pickle(
            [](const fwdpy11::GeneticValueIsFitness& g) { return g.pickle(); },
            [](py::object o) {
                std::string s = o.cast<std::string>();
                if (s.find("GeneticValueIsFitness") == std::string::npos)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                return fwdpy11::GeneticValueIsFitness();
            }));
}
