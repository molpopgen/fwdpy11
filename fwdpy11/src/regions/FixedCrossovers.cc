#include <fwdpy11/regions/FixedCrossovers.hpp>

namespace py = pybind11;

void
init_FixedCrossovers(py::module& m)
{
    py::class_<fwdpy11::FixedCrossovers, fwdpy11::GeneticMapUnit>(
        m, "FixedCrossovers",
        R"delim(
        Generate a fixed number of crossover breakpoints.
        
        .. versionadded:: 0.3.0
        )delim")
        .def(py::init<double, double, int>(), py::arg("beg"), py::arg("end"),
             py::arg("num_xovers"))
        .def_readonly("beg", &fwdpy11::FixedCrossovers::beg,
                      "Beginning of interval")
        .def_readonly("end", &fwdpy11::FixedCrossovers::end,
                      "End of interval.")
        .def_readonly("nxovers", &fwdpy11::FixedCrossovers::nxovers,
                      "Mean number of breakpoints")
        .def(py::pickle(
            [](const fwdpy11::FixedCrossovers& pi) { return pi.pickle(); },
            [](py::object o) {
                return fwdpy11::FixedCrossovers::unpickle(o);
            }));
}

