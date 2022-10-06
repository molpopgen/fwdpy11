#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_forward_graph(py::module& m);

void
init_demes(py::module& m)
{
    init_forward_graph(m);
}
