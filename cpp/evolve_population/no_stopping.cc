#include <pybind11/pybind11.h>
#include <fwdpy11/types/Population.hpp>

bool
no_stopping(const fwdpy11::Population &, const bool)
{
    return false;
}

void init_no_stopping(pybind11::module & m)
{
    m.def("_no_stopping", &no_stopping);
}

