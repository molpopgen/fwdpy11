#include "fwdpy11/regions/RecombinationRegions.hpp"
#include <pybind11/pybind11.h>
#include <core/genetic_maps/regions.hpp>

namespace py = pybind11;

void
init_GeneticMapUnit(py::module& m)
{
    py::class_<fwdpy11::PoissonCrossoverGenerator>(m, "PoissonCrossoverGenerator");
    py::class_<fwdpy11::NonPoissonCrossoverGenerator>(m, "NonPoissonCrossoverGenerator");

    py::class_<fwdpy11_core::PoissonInterval, fwdpy11::PoissonCrossoverGenerator>(
        m, "_ll_PoissonInterval")
        .def(py::init<double, double, double, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("mean"), py::arg("discrete") = true);

    py::class_<fwdpy11_core::PoissonPoint, fwdpy11::PoissonCrossoverGenerator>(
        m, "_ll_PoissonPoint")
        .def(py::init<double, double, bool>(), py::kw_only(), py::arg("position"),
             py::arg("mean"), py::arg("discrete") = true);

    py::class_<fwdpy11_core::BinomialInterval, fwdpy11::NonPoissonCrossoverGenerator>(
        m, "_ll_BinomialInterval")
        .def(py::init<double, double, double, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("probability"), py::arg("discrete") = true);

    py::class_<fwdpy11_core::BinomialPoint, fwdpy11::NonPoissonCrossoverGenerator>(
        m, "_ll_BinomialPoint")
        .def(py::init<double, double, bool>(), py::kw_only(), py::arg("position"),
             py::arg("probability"), py::arg("discrete") = true);

    py::class_<fwdpy11_core::FixedCrossovers, fwdpy11::NonPoissonCrossoverGenerator>(
        m, "_ll_FixedCrossovers")
        .def(py::init<double, double, int, bool>(), py::kw_only(), py::arg("beg"),
             py::arg("end"), py::arg("num_xovers"), py::arg("discrete") = true);
}
