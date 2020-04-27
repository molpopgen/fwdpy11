#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_data_matrix_creation_functions(py::module &);

void
initialize_fwdpp_functions(py::module &m)
{
    init_data_matrix_creation_functions(m);
}
