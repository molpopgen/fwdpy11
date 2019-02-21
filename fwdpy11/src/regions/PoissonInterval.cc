#include <fwdpy11/regions/PoissonInterval.hpp>

namespace py = pybind11;

void
init_PoissonInterval(py::module& m)
{
    py::class_<fwdpy11::PoissonInterval, fwdpy11::GeneticMapUnit>(
        m, "PoissonInterval",
        R"delim(
        Generate poisson number of crossover breakpoints.
        
        .. versionadded:: 0.3.0
        )delim")
        .def(py::init<double, double, double>(), py::arg("beg"),
             py::arg("end"), py::arg("mean"))
        .def_readonly("beg", &fwdpy11::PoissonInterval::beg,
                      "Beginning of interval")
        .def_readonly("end", &fwdpy11::PoissonInterval::end,
                      "End of interval.")
        .def_readonly("mean", &fwdpy11::PoissonInterval::mean,
                      "Mean number of breakpoints")
        .def(py::pickle(
            [](const fwdpy11::PoissonInterval& pi) { return pi.pickle(); },
            [](py::object o) {
                return fwdpy11::PoissonInterval::unpickle(o);
            }));
}

