#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GeneticMapUnit.hpp>

namespace py = pybind11;

void
init_GeneticMapUnit(py::module& m)
{
    py::class_<fwdpy11::GeneticMapUnit>(m, "GeneticMapUnit");

    py::class_<fwdpy11::PoissonInterval, fwdpy11::GeneticMapUnit>(m,
                                                                  "_ll_PoissonInterval")
        .def(py::init<double, double, double, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("mean"), py::arg("discrete") = true);

    py::class_<fwdpy11::PoissonPoint, fwdpy11::GeneticMapUnit>(m, "_ll_PoissonPoint")
        .def(py::init<double, double, bool>(), py::kw_only(), py::arg("position"),
             py::arg("mean"), py::arg("discrete") = true);

    py::class_<fwdpy11::BinomialInterval, fwdpy11::GeneticMapUnit>(
        m, "_ll_BinomialInterval")
        .def(py::init<double, double, double, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("probability"), py::arg("discrete") = true);

    py::class_<fwdpy11::BinomialPoint, fwdpy11::GeneticMapUnit>(m, "_ll_BinomialPoint")
        .def(py::init<double, double, bool>(), py::kw_only(), py::arg("position"),
             py::arg("probability"), py::arg("discrete") = true);

    py::class_<fwdpy11::FixedCrossovers, fwdpy11::GeneticMapUnit>(m,
                                                                  "_ll_FixedCrossovers")
        .def(py::init<double, double, int, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("num_xovers"), py::arg("discrete") = true);
}
