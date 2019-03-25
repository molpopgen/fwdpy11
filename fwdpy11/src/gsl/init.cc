#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_gsl_random(py::module&);

void init_GSL(py::module & m)
{
    init_gsl_random(m);
}
