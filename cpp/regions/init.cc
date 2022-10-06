#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Region(py::module &);
void init_Sregion(py::module &);
void init_GammaS(py::module &);
void init_ConstantS(py::module &);
void init_ExpS(py::module &);
void init_UniformS(py::module &);
void init_GaussianS(py::module &);
void init_MutationRegions(py::module &);
void init_RecombinationRegions(py::module &);
void init_MultivariateGaussianEffects(py::module &);
void init_mvDES(py::module &);
void init_LogNormalS(py::module &);
void init_DiscreteDESD(py::module&);

void
initialize_regions(py::module &m)
{
    init_Region(m);
    init_Sregion(m);
    init_GammaS(m);
    init_ConstantS(m);
    init_ExpS(m);
    init_UniformS(m);
    init_GaussianS(m);
    init_MutationRegions(m);
    init_RecombinationRegions(m);
    init_MultivariateGaussianEffects(m);
    init_mvDES(m);
    init_LogNormalS(m);
    init_DiscreteDESD(m);
}
