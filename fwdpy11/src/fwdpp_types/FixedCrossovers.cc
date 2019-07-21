#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/fixed_number_crossovers.hpp>

namespace py = pybind11;

void
init_FixedCrossovers(py::module& m)
{
    py::class_<fwdpp::fixed_number_crossovers, fwdpp::genetic_map_unit>(
        m, "FixedCrossovers",
        R"delim(
        Generate a fixed number of crossover breakpoints.
        
        .. versionadded:: 0.3.0

        .. versionchanged:: 0.5.0

            Refactored back-end to be based on fwdpp types
        )delim")
        .def(py::init<double, double, int>(), py::arg("beg"), py::arg("end"),
             py::arg("num_xovers"))
        .def_readonly("beg", &fwdpp::fixed_number_crossovers::beg,
                      "Beginning of interval")
        .def_readonly("end", &fwdpp::fixed_number_crossovers::end,
                      "End of interval.")
        .def_readonly("nxovers", &fwdpp::fixed_number_crossovers::nxovers,
                      "Mean number of breakpoints")
        .def(py::pickle(
            [](const fwdpp::fixed_number_crossovers& self) {
                return py::make_tuple(self.beg, self.end, self.nxovers);
            },
            [](py::tuple t) {
                return std::unique_ptr<fwdpp::fixed_number_crossovers>(
                    new fwdpp::fixed_number_crossovers(t[0].cast<double>(),
                                                       t[1].cast<double>(),
                                                       t[2].cast<int>()));
            }));
}

