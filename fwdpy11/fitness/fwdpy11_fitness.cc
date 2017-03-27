#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fwdpy11_fitness)
{
    py::module m("fwdpy11_fitness", "Standard fitness models.");

	py::class_<fwdpy11::singlepop_fitness>(m,"SpopFitness");

    py::class_<fwdpy11::singlepop_mult_wrapper,fwdpy11::singlepop_fitness>(m, "SpopMult")
        .def(py::init<double>());

    py::class_<fwdpy11::singlepop_additive_wrapper,fwdpy11::singlepop_fitness>(m, "SpopAdditive")
        .def(py::init<double>());

    return m.ptr();
}
