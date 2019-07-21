#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/poisson_interval.hpp>

namespace py = pybind11;

void
init_PoissonInterval(py::module& m)
{
    py::class_<fwdpp::poisson_interval, fwdpp::genetic_map_unit>(
        m, "PoissonInterval",
        R"delim(
        Generate poisson number of crossover breakpoints.
        
        .. versionadded:: 0.3.0

        .. versionchanged:: 0.5.0

            Refactored back-end to be based on fwdpp types
        )delim")
        .def(py::init<double, double, double>(), py::arg("beg"),
             py::arg("end"), py::arg("mean"))
        .def_readonly("beg", &fwdpp::poisson_interval::beg,
                      "Beginning of interval")
        .def_readonly("end", &fwdpp::poisson_interval::end, "End of interval.")
        .def_readonly("mean", &fwdpp::poisson_interval::mean,
                      "Mean number of breakpoints")
        .def(py::pickle(
            [](const fwdpp::poisson_interval& self) {
                return py::make_tuple(self.beg, self.end, self.mean);
            },
            [](py::tuple t) {
                return std::unique_ptr<fwdpp::poisson_interval>(
                    new fwdpp::poisson_interval(t[0].cast<double>(),
                                                t[1].cast<double>(),
                                                t[2].cast<double>()));
            }));
}

