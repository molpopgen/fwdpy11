#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Region(py::module &);
void init_Sregion(py::module &);
void init_GammaS(py::module &);
void init_ConstantS(py::module &);
void init_ExpS(py::module &);
void init_UniformS(py::module &);
void init_GaussianS(py::module &);

void initialize_regions(py::module & m)
{
    init_Region(m);
    init_Sregion(m);
    init_GammaS(m);
    init_ConstantS(m);
    init_ExpS(m);
    init_UniformS(m);
    init_GaussianS(m);
}

PYBIND11_MODULE(_fwdpy11, m)
{
    initialize_regions(m);
}
