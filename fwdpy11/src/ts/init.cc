#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_tree_iterator(py::module&);

void
init_ts(py::module& m)
{
    init_tree_iterator(m);
}
