#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include "types.hpp"

namespace py = pybind11;

using singlepop_base = fwdpy::singlepop_t::popbase;

PYBIND11_PLUGIN(fwdpy11_types) {
    py::module m("fwdpy11_types", "example extending");

    py::class_<fwdpy::GSLrng_t>(
        m, "GSLrng", "Random number generator based on a mersenne twister.")
        .def(py::init<unsigned>(),
             "Constructor takes unsigned integer as a seed");

    py::class_<fwdpy::diploid_t>(
        m, "SingleLocusDiploid",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_readonly("first", &fwdpy::diploid_t::first)
        .def_readonly("second", &fwdpy::diploid_t::second)
        .def_readonly("w", &fwdpy::diploid_t::w)
        .def_readonly("g", &fwdpy::diploid_t::g)
        .def_readonly("e", &fwdpy::diploid_t::e);

    pybind11::class_<singlepop_base>(m, "SinglepopBase");

    // Expose the type based on fwdpp's "sugar" layer
    pybind11::class_<fwdpy::singlepop_t, singlepop_base>(m, "Spop")
        .def(pybind11::init<unsigned>())
        .def("clear", &fwdpy::singlepop_t::clear)
        .def_readonly("mutations", &fwdpy::singlepop_t::mutations)
        .def_readonly("mcounts", &fwdpy::singlepop_t::mcounts)
        .def_readonly("fixations", &fwdpy::singlepop_t::fixations)
        .def_readonly("diploids", &fwdpy::singlepop_t::diploids)
        .def_readonly("gametes", &fwdpy::singlepop_t::gametes)
        .def_readonly("generation", &fwdpy::singlepop_t::generation)
        .def_readonly("N", &fwdpy::singlepop_t::N)
        .def("__getstate__",
             [](const fwdpy::singlepop_t& pop) {
                 return py::make_tuple(py::bytes(pop.serialize()));
             })
        .def("__setstate__", [](fwdpy::singlepop_t& p, py::tuple s) {
            new (&p) fwdpy::singlepop_t(s[0].cast<std::string>());
        });

    return m.ptr();
}
