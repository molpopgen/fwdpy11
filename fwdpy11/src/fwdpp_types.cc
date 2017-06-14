//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>
#include <fwdpy11/opaque/opaque_types.hpp>
namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::uint_t>);
PYBIND11_MAKE_OPAQUE(
    std::vector<double>); // for generalmut_vec::s and generalmut_vec::h

PYBIND11_PLUGIN(fwdpp_types)
{
    py::module m("fwdpp_types", "Wrap C++ types from fwdpp.");

    // low-level types

    py::class_<KTfwd::mutation_base>(m, "MutationBase",
                                     R"delim(
Base class for mutations.
)delim")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readwrite("pos", &KTfwd::mutation_base::pos, "Position (float).")
        .def_readwrite("neutral", &KTfwd::mutation_base::neutral, "Boolean")
        .def_readwrite("xtra", &KTfwd::mutation_base::xtra,
                       "16-bits worth of extra data.");

    py::class_<KTfwd::gamete>(m, "Gamete", R"delim(
    A gamete.  This object represents a haplotype
    in a contiguous genomic region.
)delim")
        .def_readonly("n", &KTfwd::gamete::n,
                      "Number of occurrences in the population.  . This has "
                      "little meaning beyond book-keeping used by the C++ "
                      "back-end. (read-only)")
        .def_readonly("mutations", &KTfwd::gamete::mutations,
                      "List of keys to neutral mutations. Contains unsigned "
                      "32-bit integers corresponding to mutations in the "
                      "population. (read-only)")
        .def_readonly("smutations", &KTfwd::gamete::smutations,
                      "List of keys to selected mutations. Contains unsigned "
                      "32-bit integers corresponding to mutations in the "
                      "population. (read-only)")
        .def("as_dict",
             // This lambda shows that
             // we can return dicts
             // with a mix of scalars + containers
             [](const KTfwd::gamete &g) noexcept {
                 using obj = pybind11::object;
                 pybind11::dict rv;
                 rv[obj(pybind11::cast("n"))] = obj(pybind11::cast(g.n));
                 rv[obj(pybind11::cast("mutations"))]
                     = obj(pybind11::cast(g.mutations));
                 rv[obj(pybind11::cast("smutations"))]
                     = obj(pybind11::cast(g.smutations));
                 return rv;
             },
             "Return dictionary representaton of the gamete.")
        .def("__getstate__",
             [](const KTfwd::gamete &g) {
                 return py::make_tuple(g.n, g.mutations, g.smutations);
             })
        .def("__setstate__", [](KTfwd::gamete &g, py::tuple t) {
            new (&g) KTfwd::gamete(t[0].cast<KTfwd::uint_t>(),
                                   t[1].cast<std::vector<KTfwd::uint_t>>(),
                                   t[2].cast<std::vector<KTfwd::uint_t>>());
        });

    // Sugar types
    py::class_<KTfwd::popgenmut, KTfwd::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init<double, double, double, unsigned, std::uint16_t>())
        .def_readwrite(
            "g", &KTfwd::popgenmut::g,
            "Generation when mutation arose (origination time). (read-only)")
        .def_readwrite("s", &KTfwd::popgenmut::s,
                       "Selection coefficient/effect size. (read-only)")
        .def_readwrite("h", &KTfwd::popgenmut::h,
                       "Dominance/effect in heterozygotes. (read-only)")
        .def("__getstate__",
             [](const KTfwd::popgenmut &m) {
                 return py::make_tuple(m.pos, m.s, m.h, m.g, m.xtra);
             })
        .def("__setstate__",
             [](KTfwd::popgenmut &m, py::tuple p) {
                 new (&m) KTfwd::popgenmut(
                     p[0].cast<double>(), p[1].cast<double>(),
                     p[2].cast<double>(), p[3].cast<unsigned>(),
                     p[4].cast<std::uint16_t>());
             })
        .def("__str__", [](const KTfwd::popgenmut &m) {
            return "Mutation[" + std::to_string(m.pos) + ","
                   + std::to_string(m.s) + "," + std::to_string(m.h) + ","
                   + std::to_string(m.g) + "]";
        });

    py::bind_vector<std::vector<double>>(
        m, "VectorDouble", "Vector of 64-bit floats.", py::buffer_protocol());
    py::bind_vector<std::vector<KTfwd::generalmut_vec>>(
        m, "VectorGeneralMutVec",
        "A list of :class:`fwdpy11.fwdpp_types.GeneralMutVec`.");

    py::class_<KTfwd::generalmut_vec, KTfwd::mutation_base>(
        m, "GeneralMutVec",
        "Mutation type with vector of effect size and dominance terms.")
        .def_readonly("s", &KTfwd::generalmut_vec::s,
                      "List of selection coefficients/effect sizes.")
        .def_readonly("h", &KTfwd::generalmut_vec::h,
                      "List of dominance terms.")
        .def_readonly("g", &KTfwd::generalmut_vec::g,
                      "Generation when mutation arose.");

    return m.ptr();
}
