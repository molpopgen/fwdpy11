#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <type_traits>
#include "types.hpp"

namespace py = pybind11;

using fwdpp_popgenmut_base = fwdpy::singlepop_t::popbase_t;
using singlepop_sugar_base = fwdpy::singlepop_t::base;

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
        .def_readonly("first", &fwdpy::diploid_t::first, "Key to first gamete.")
        .def_readonly("second", &fwdpy::diploid_t::second,
                      "Key to second gamete.")
        .def_readonly("w", &fwdpy::diploid_t::w, "Fitness.")
        .def_readonly("g", &fwdpy::diploid_t::g, "Genetic value.")
        .def_readonly("e", &fwdpy::diploid_t::e,
                      "Random/environmental effects.")
        .def("__getstate__",
             [](const fwdpy::diploid_t& d) {
                 return py::make_tuple(d.first, d.second, d.w, d.g, d.e);
             })
        .def("__setstate__", [](fwdpy::diploid_t& d, py::tuple t) {
            new (&d) fwdpy::diploid_t(t[0].cast<std::size_t>(),
                                      t[1].cast<std::size_t>());
            d.w = t[2].cast<double>();
            d.g = t[3].cast<double>();
            d.e = t[4].cast<double>();
        });

    py::class_<fwdpp_popgenmut_base>(m, "MutationPoptypeCommonBase")
        .def_readonly("mutations", &fwdpp_popgenmut_base::mutations,
                      "Container of :class:`fwdpy11.fwdpp_types.Mutation`")
        .def_readonly("mcounts", &fwdpp_popgenmut_base::mcounts)
        .def_readonly("fixations", &fwdpp_popgenmut_base::fixations)
        .def_readonly("gametes", &fwdpp_popgenmut_base::gametes);

    py::class_<singlepop_sugar_base, fwdpp_popgenmut_base>(m, "SinglepopBase")
        .def_readonly("diploids", &singlepop_sugar_base::diploids);

    // Expose the type based on fwdpp's "sugar" layer
    py::class_<fwdpy::singlepop_t, singlepop_sugar_base>(
        m, "Spop", "Population object representing a single deme.")
        .def(py::init<unsigned>(),
             "Construct with an unsigned integer representing the initial "
             "population size.")
        .def("clear", &fwdpy::singlepop_t::clear, "Clears all population data.")
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
