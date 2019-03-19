#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_change_effect_size(py::module &);

void initialize_fwdpy11_functions(py::module & m)
{
    init_change_effect_size(m);
}
