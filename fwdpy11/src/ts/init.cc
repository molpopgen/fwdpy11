#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_tree_iterator(py::module&);
void init_variant_iterator(py::module&);
void init_count_mutations(py::module&);
void init_simplify_functions(py::module&);
void init_data_matrix_from_tables(py::module&);
void init_infinite_sites(py::module&);
void
init_DataMatrixIterator(py::module& m);

void
init_ts(py::module& m)
{
    init_tree_iterator(m);
    init_variant_iterator(m);
    init_count_mutations(m);
    init_simplify_functions(m);
    init_data_matrix_from_tables(m);
    init_infinite_sites(m);
    init_DataMatrixIterator(m);
}
