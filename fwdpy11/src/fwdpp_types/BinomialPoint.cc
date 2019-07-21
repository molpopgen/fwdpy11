#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/binomial_point.hpp>

namespace py = pybind11;

void
init_BinomialPoint(py::module& m)
{
    py::class_<fwdpp::binomial_point, fwdpp::genetic_map_unit>(m,
                                                               "BinomialPoint",
                                                               R"delim(
        Generate a crossover breakpoint at a fixed position with a
        fixed probability.  This class represents genetic distance
        as centiMorgans/100.

        .. versionadded:: 0.3.0

        .. versionchanged:: 0.5.0

            Refactored back-end to be based on fwdpp types
        )delim")
        .def(py::init<double, double>(), py::arg("position"),
             py::arg("probability"))
        .def_readonly("position", &fwdpp::binomial_point::position, "Position")
        .def_readonly("probability", &fwdpp::binomial_point::prob,
                      "Probability")
        .def(py::pickle(
            [](const fwdpp::binomial_point& self) {
                return py::make_tuple(self.position, self.prob);
            },
            [](py::tuple t) {
                return std::unique_ptr<fwdpp::binomial_point>(
                    new fwdpp::binomial_point(t[0].cast<double>(),
                                              t[1].cast<double>()));
            }));
}

