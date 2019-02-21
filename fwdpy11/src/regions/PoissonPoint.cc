#include <fwdpy11/regions/PoissonPoint.hpp>

namespace py = pybind11;

void
init_PoissonPoint(py::module& m)
{
    py::class_<fwdpy11::PoissonPoint, fwdpy11::GeneticMapUnit>(m,
                                                               "PoissonPoint")
        .def(py::init<double, double>(), py::arg("position"), py::arg("rate"))
        .def_readonly("position", &fwdpy11::PoissonPoint::position)
        .def_readonly("rate", &fwdpy11::PoissonPoint::rate)
        .def(py::pickle(
            [](const fwdpy11::PoissonPoint& b) { return b.pickle(); },
            [](py::object o) { return fwdpy11::PoissonPoint::unpickle(o); }));
}

