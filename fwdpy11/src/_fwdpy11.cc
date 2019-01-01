#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Region(py::module &);
void init_Sregion(py::module &);
void init_GammaS(py::module &);
void init_ConstantS(py::module &);
void init_ExpS(py::module &);

PYBIND11_MODULE(_fwdpy11, m)
{
    //Regions
    init_Region(m);
    init_Sregion(m);
    init_GammaS(m);
    init_ConstantS(m);
    init_ExpS(m);
}
