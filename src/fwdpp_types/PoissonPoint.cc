#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/poisson_point.hpp>

namespace py = pybind11;
using namespace py::literals;

void
init_PoissonPoint(py::module& m)
{
    py::class_<fwdpp::poisson_point, fwdpp::genetic_map_unit>(m, "PoissonPoint",
                                                              R"delim(
        Generate a recombination breakpoint at a fixed position if the
        number of crossover events is odd.

        .. versionadded:: 0.3.0

        .. versionchanged:: 0.5.0

            Refactored back-end to be based on fwdpp types

        .. versionchanged:: 0.7.1

            Add __repr__
        )delim")
        .def(py::init<double, double>(), py::arg("position"), py::arg("mean"))
        .def_readonly("position", &fwdpp::poisson_point::position, "Position")
        .def_readonly("mean", &fwdpp::poisson_point::mean, "Mean of Poisson process")
        .def("__repr__",
             [](const fwdpp::poisson_point& self) {
                 return "PoissonPoint(position={}, mean={})"_s.format(self.position,
                                                                      self.mean);
             })
        .def(py::pickle(
            [](const fwdpp::poisson_point& self) {
                return py::make_tuple(self.position, self.mean);
            },
            [](py::tuple t) {
                return std::unique_ptr<fwdpp::poisson_point>(
                    new fwdpp::poisson_point(t[0].cast<double>(), t[1].cast<double>()));
            }));
}

