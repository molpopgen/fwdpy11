#include <pybind11/pybind11.h>

namespace py = pybind11;

void initialize_regions(py::module &);

PYBIND11_MODULE(_fwdpy11, m)
{
    initialize_regions(m);
}
