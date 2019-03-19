#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

void init_no_stopping(py::module &);
void init_diploid_evolution(py::module &);

PYBIND11_MODULE(_tsevolution, m)
{
    m.doc() = "Evolution under a Wright-Fisher model using tree sequences.";

    init_no_stopping(m);
    init_diploid_evolution(m);
}

