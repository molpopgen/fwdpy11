#include <pybind11/pybind11.h>
#include <fwdpp/genetic_map/genetic_map_unit.hpp>

namespace py = pybind11;

void
init_GeneticMapUnit(py::module& m)
{
    py::class_<fwdpp::genetic_map_unit>(m, "GeneticMapUnit",
                                        R"delim(
                                        ABC for components of genetic maps

                                        .. versionadded:: 0.3.0
                                        
                                        .. versionchanged:: 0.5.0

                                            Refactored back-end to be based on fwdpp types
                                        )delim");
}
