#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fwdpy11_fitness)
{
    py::module m("fwdpy11_fitness", "Standard fitness models.");

    py::class_<fwdpy11::singlepop_mult_wrapper>(m, "SpopMult")
        .def(py::init<double>())
        .def_property_readonly(
            "callback",
            [](const fwdpy11::singlepop_mult_wrapper &m) -> fwdpy11::singlepop_fitness {
                return m.callback();
            });

    py::class_<fwdpy11::singlepop_additive_wrapper>(m, "SpopAdditive")
        .def(py::init<double>())
        .def_property_readonly(
            "callback",
            [](const fwdpy11::singlepop_additive_wrapper &m) -> fwdpy11::singlepop_fitness {
                return m.callback();
            });

    return m.ptr();
}
