#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

void init_no_stopping(py::module &);
void init_evolve_with_tree_sequences(py::module &);

void
init_evolution_functions(py::module &m)
{
    init_no_stopping(m);
    init_evolve_with_tree_sequences(m);
}

