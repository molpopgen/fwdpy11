#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/binomial_interval.hpp>

namespace py = pybind11;
using namespace py::literals;

void
init_BinomialInterval(py::module& m)
{
    py::class_<fwdpp::binomial_interval, fwdpp::genetic_map_unit>(m, "BinomialInterval",
                                                                  R"delim(
        Generate exactly one crossover with a given probability
        
        .. versionadded:: 0.5.2

        .. versionchanged:: 0.7.1

            Add __repr__
        )delim")
        .def(py::init<double, double, double>(), py::arg("beg"), py::arg("end"),
             py::arg("probability"))
        .def_readonly("beg", &fwdpp::binomial_interval::beg, "Beginning of interval")
        .def_readonly("end", &fwdpp::binomial_interval::end, "End of interval.")
        .def_readonly("probability", &fwdpp::binomial_interval::prob,
                      "Probability of a crossover")
        .def("__repr__",
             [](const fwdpp::binomial_interval& self) {
                 return "BinomialInterval(beg={}, end={}, probability={})"_s.format(
                     self.beg, self.end, self.prob);
             })
        .def(py::pickle(
            [](const fwdpp::binomial_interval& self) {
                return py::make_tuple(self.beg, self.end, self.prob);
            },
            [](py::tuple t) {
                return std::unique_ptr<fwdpp::binomial_interval>(
                    new fwdpp::binomial_interval(
                        t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()));
            }));
}

