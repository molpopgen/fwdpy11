#include <fwdpy11/regions/PoissonPoint.hpp>

namespace py = pybind11;

void
init_PoissonPoint(py::module& m)
{
    py::class_<fwdpy11::PoissonPoint, fwdpy11::GeneticMapUnit>(m,
                                                               "PoissonPoint",
                                                               R"delim(
        Generate a recombination breakpoint at a fixed position if the
        number of crossover events is odd.

        .. versionadded:: 0.3.0
        )delim")
        .def(py::init<double, double>(), py::arg("position"), py::arg("mean"))
        .def_readonly("position", &fwdpy11::PoissonPoint::position, "Position")
        .def_readonly("mean", &fwdpy11::PoissonPoint::mean,
                      "Mean of Poisson process")
        .def(py::pickle(
            [](const fwdpy11::PoissonPoint& b) { return b.pickle(); },
            [](py::object o) { return fwdpy11::PoissonPoint::unpickle(o); }));
}

