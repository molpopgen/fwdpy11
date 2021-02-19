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

    py::class_<fwdpy11::FixedDominance, fwdpy11::MutationDominance>(m,
                                                                    "_ll_FixedDominance")
        .def(py::init<double>(), py::arg("h"));
    py::class_<fwdpy11::ExponentialDominance, fwdpy11::MutationDominance>(
        m, "_ll_ExponentialDominance")
        .def(py::init<double>(), py::arg("m"));
    py::class_<fwdpy11::UniformDominance, fwdpy11::MutationDominance>(
        m, "_ll_UniformDominance")
        .def(py::init<double, double>(), py::arg("lo"), py::arg("hi"));
    py::class_<fwdpy11::LargeEffectExponentiallyRecessive, fwdpy11::MutationDominance>(
        m, "_ll_LargeEffectExponentiallyRecessive")
        .def(py::init<double, double>(), py::arg("k"), py::arg("scaling"));
}
