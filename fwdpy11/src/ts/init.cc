#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_tree_iterator(py::module&);
void init_variant_iterator(py::module&);
void init_count_mutations(py::module&);

void
init_ts(py::module& m)
{
    init_tree_iterator(m);
    init_variant_iterator(m);
    init_count_mutations(m);
}
