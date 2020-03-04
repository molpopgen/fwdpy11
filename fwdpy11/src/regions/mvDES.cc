#include <pybind11/pybind11.h>
#include <fwdpy11/regions/mvDES.hpp>

namespace py = pybind11;

void
init_mvDES(py::module &m)
{
    py::class_<fwdpy11::mvDES, fwdpy11::Sregion>(m, "mvDES");
}
