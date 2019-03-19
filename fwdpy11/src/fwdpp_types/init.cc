#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_mutation_base(py::module &);
void init_gamete(py::module &);

void
initialize_fwdpp_types(py::module &m)
{
    init_mutation_base(m);
    init_gamete(m);
}
