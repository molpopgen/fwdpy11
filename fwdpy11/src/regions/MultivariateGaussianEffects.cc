#include <pybind11/pybind11.h>
#include <fwdpy11/regions/MultivariateGaussianEffects.hpp>

namespace py = pybind11;

void
init_MultivariateGaussianEffects(py::module& m)
{
    py::class_<fwdpy11::MultivariateGaussianEffects, fwdpy11::Sregion>(
        m, "MultivariateGaussianEffects");
}
