#include <fwdpy11/regions/MutationDominance.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_MutationDominance(py::module& m)
{
    py::class_<fwdpy11::MutationDominance>(m, "MutationDominance",
                                           R"delim(
            ABC for classes determining the dominance of new mutations.

            .. versionadded:: 0.13.0
            )delim");
}
