#include <fwdpy11/regions/PoissonInterval.hpp>

namespace py = pybind11;

void
init_PoissonInterval(py::module& m)
{
    py::class_<fwdpy11::PoissonInterval, fwdpy11::GeneticMapUnit>(
        m, "PoissonInterval")
        .def(py::init<double, double, double>(), py::arg("beg"),
             py::arg("end"), py::arg("rate"))
        .def_readonly("beg", &fwdpy11::PoissonInterval::beg)
        .def_readonly("end", &fwdpy11::PoissonInterval::end)
        .def_readonly("rate", &fwdpy11::PoissonInterval::rate)
        .def(py::pickle(
            [](const fwdpy11::PoissonInterval& pi) { return pi.pickle(); },
            [](py::object o) {
                return fwdpy11::PoissonInterval::unpickle(o);
            }));
}

