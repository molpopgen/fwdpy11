#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/regions/RecombinationRegions.hpp>

namespace py = pybind11;

void
init_RecombinationRegions(py::module& m)
{
    py::class_<fwdpy11::GeneticMap>(m, "GeneticMap", "ABC for genetic maps");

    py::class_<fwdpy11::RecombinationRegions, fwdpy11::GeneticMap>(
        m, "RecombinationRegions")
        .def(py::init<double, std::vector<fwdpy11::Region>>())
        .def_readonly("weights", &fwdpy11::RecombinationRegions::weights);

    py::class_<fwdpy11::MlocusRecombinationRegions>(
        m, "MlocusRecombinationRegions")
        .def(py::init<>())
        .def("append", &fwdpy11::MlocusRecombinationRegions::append);
}
