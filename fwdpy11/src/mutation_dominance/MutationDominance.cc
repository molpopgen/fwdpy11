#include <fwdpy11/mutation_dominance/MutationDominance.hpp>
#include <fwdpy11/mutation_dominance/ExponentialDominance.hpp>
#include <fwdpy11/mutation_dominance/UniformDominance.hpp>
#include <fwdpy11/mutation_dominance/LargeEffectExponentiallyRecessive.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_MutationDominance(py::module& m)
{
    py::class_<fwdpy11::MutationDominance>(m, "MutationDominance",
                                           R"delim(
            ABC for classes determining the dominance of new mutations.

            .. versionadded:: 0.13.0
            )delim");

    py::class_<fwdpy11::FixedDominance>(m, "_ll_FixedDominance").def(py::init<double>());
    py::class_<fwdpy11::ExponentialDominance>(m, "_ll_ExponentialDominance")
        .def(py::init<double>());
    py::class_<fwdpy11::UniformDominance>(m, "_ll_UniformDominance")
        .def(py::init<double, double>());
    py::class_<fwdpy11::LargeEffectExponentiallyRecessive>(m, "_ll_LargeEffectExponentiallyRecessive")
        .def(py::init<double>());
}