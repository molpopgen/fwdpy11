#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

void init_no_stopping(py::module &);
void init_evolve_with_tree_sequences(py::module &);
void init_evolve_without_tree_sequences(py::module & m);

PYBIND11_MODULE(_evolve_population, m)
{
    m.doc() = "Evolution under a Wright-Fisher model using tree sequences.";

    init_no_stopping(m);
    init_evolve_with_tree_sequences(m);
    init_evolve_without_tree_sequences(m);
}

