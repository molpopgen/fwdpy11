#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fitness)
{
    py::module m("fitness", "Fitness models.");

    py::class_<fwdpy11::singlepop_fitness>(m, "SpopFitness");
    py::class_<fwdpy11::singlepop_fitness_qtrait>(m, "SpopFitnessQtrait");

    py::class_<fwdpy11::singlepop_mult_wrapper, fwdpy11::singlepop_fitness>(
        m, "SpopMult", R"delim(
        Multiplicative fitness for single-deme simulations.
        At a single mutation, fitness is 0, 1+sh, 1+scaling*s
        for genotypes AA, Aa, and aa, respectively. The scaling
        parameter is a constructor argument.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.SpopMult(1.0)
                       )delim")
        .def(py::init<double>(), py::arg("scaling"));

    py::class_<fwdpy11::singlepop_additive_wrapper,
               fwdpy11::singlepop_fitness>(
        m, "SpopAdditive", R"delim(
        Additive fitness for single-deme simulations.
        Fitness is max(0,1 + :math:`\sum_{i} x_i`), 
        where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
        for AA, Aa, and aa, respectively.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.SpopAdditive(2.0)
        )delim")
        .def(py::init<double>(), py::arg("scaling"));

    return m.ptr();
}
