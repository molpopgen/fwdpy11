#include <fwdpy11/mutation_dominance/MutationDominance.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_MutationDominance(py::module& m)
{
    py::class_<fwdpy11::MutationDominance>(m, "MutationDominance",
                                           R"delim(
            ABC for classes determining the dominance of new mutations.

            .. versionadded:: 0.13.0
            )delim")
        .def(py::init<fwdpy11::MutationDominance>());

    m.def("_fixed_dominance", &fwdpy11::fixed_dominance);
    m.def("_uniform_dominance", &fwdpy11::uniform_dominance);
    m.def("_exponential_dominance", &fwdpy11::exponential_dominance);
    m.def("_large_effect_exponentially_recessive",
          &fwdpy11::large_effect_exponentially_recessive);
}
