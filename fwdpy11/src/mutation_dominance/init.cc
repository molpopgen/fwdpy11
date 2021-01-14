#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_MutationDominance(py::module &);

void initialize_mutation_dominance(py::module &m)
{
    init_MutationDominance(m);
}
