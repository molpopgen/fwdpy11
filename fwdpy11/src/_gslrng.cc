#include <pybind11/pybind11.h>
#include <fwdpy11/rng.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_gslrng, m)
{
    py::class_<fwdpy11::GSLrng_t>(m, "GSLrng", "Random number generator based "
                                               "on GNU Scientific Library "
                                               "mersenne twister.")
        .def(py::init<unsigned>(),
             "Constructor takes unsigned integer as a seed");
}
