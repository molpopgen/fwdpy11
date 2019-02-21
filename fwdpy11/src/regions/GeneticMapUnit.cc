#include <fwdpy11/regions/GeneticMapUnit.hpp>

namespace py = pybind11;

void
init_GeneticMapUnit(py::module& m)
{
    py::class_<fwdpy11::GeneticMapUnit>(m, "GeneticMapUnit",
                                        "ABC for components of genetic maps");
}
