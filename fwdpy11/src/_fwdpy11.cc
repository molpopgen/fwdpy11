#include <pybind11/pybind11.h>

namespace py = pybind11;

void initialize_fwdpp_types(py::module &);
void initialize_fwdpp_functions(py::module &);
void initialize_fwdpy11_types(py::module & m);
void initialize_regions(py::module &);
void init_GSL(py::module &);

PYBIND11_MODULE(_fwdpy11, m)
{
    initialize_fwdpp_types(m);
    initialize_fwdpp_functions(m);
    initialize_fwdpy11_types(m);
    initialize_regions(m);
    init_GSL(m);
}
