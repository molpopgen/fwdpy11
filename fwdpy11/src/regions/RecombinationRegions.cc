#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/regions/RecombinationRegions.hpp>

namespace py = pybind11;

void
init_RecombinationRegions(py::module& m)
{
    py::class_<fwdpy11::RecombinationRegions>(m, "RecombinationRegions")
        .def(py::init<double, std::vector<fwdpy11::Region>>())
        .def_readonly("weights", &fwdpy11::RecombinationRegions::weights);
}
