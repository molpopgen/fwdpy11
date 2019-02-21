#include <fwdpy11/regions/BinomialPoint.hpp>

namespace py = pybind11;

void
init_BinomialPoint(py::module& m)
{
    py::class_<fwdpy11::BinomialPoint, fwdpy11::GeneticMapUnit>(
        m, "BinomialPoint",
        R"delim(
        Generate a crossover breakpoint at a fixed position with a
        fixed probability.  This class represents genetic distance
        as centiMorgans/100.

        .. versionadded:: 0.3.0
        )delim")
        .def(py::init<double, double>(), py::arg("position"),
             py::arg("probability"))
        .def_readonly("position", &fwdpy11::BinomialPoint::position,
                      "Position")
        .def_readonly("probability", &fwdpy11::BinomialPoint::prob,
                      "Probability")
        .def(py::pickle(
            [](const fwdpy11::BinomialPoint& b) { return b.pickle(); },
            [](py::object o) { return fwdpy11::BinomialPoint::unpickle(o); }));
}

