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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::uint_t>);

PYBIND11_MODULE(fwdpp_types, m)
{
    m.doc() = "Wrap C++ types from fwdpp.";

    // low-level types

    py::class_<fwdpp::mutation_base>(m, "MutationBase",
                                     R"delim(
                                        Base class for mutations.
                                     )delim")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readonly("pos", &fwdpp::mutation_base::pos, "Position (float).")
        .def_readonly("neutral", &fwdpp::mutation_base::neutral, "Boolean")
        .def_readwrite("label", &fwdpp::mutation_base::xtra,
                       "A 16-bit unsigned integer that can be used for adding "
                       "\"meta-data\" to mutations");

    py::bind_vector<std::vector<fwdpp::uint_t>>(
        m, "VecUint32", "Vector of unsigned 32-bit integers.",
        py::buffer_protocol())
        .def(py::pickle(
            [](const std::vector<fwdpp::uint_t>& v) {
                py::list rv;
                for (auto&& i : v)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpp::uint_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpp::uint_t>());
                    }
                return rv;
            }));

    py::class_<fwdpp::gamete>(m, "Gamete", R"delim(
    A gamete.  This object represents a haplotype
    in a contiguous genomic region.
)delim")
        .def(py::init<fwdpp::gamete::constructor_tuple>(),
             R"delim(
                Construct gamete from tuple.
                
                The tuple must be (n, mutations, smutations)

                .. testcode::

                    import fwdpy11
                    # Note the cast that is needed: 
                    g = fwdpy11.Gamete((1,
                                        fwdpy11.VecUint32([2]),
                                        fwdpy11.VecUint32([0])))
                    print(g.n)
                    print(list(g.mutations))
                    print(list(g.smutations))

                .. testoutput::

                    1
                    [2]
                    [0]
                )delim")
        .def_readonly("n", &fwdpp::gamete::n,
                      "Number of occurrences in the population. This has "
                      "little meaning beyond book-keeping used by the C++ "
                      "back-end. (read-only)")
        .def_readonly("mutations", &fwdpp::gamete::mutations,
                      "List of keys to neutral mutations. Contains unsigned "
                      "32-bit integers corresponding to mutations in the "
                      "population. (read-only)")
        .def_readonly("smutations", &fwdpp::gamete::smutations,
                      "List of keys to selected mutations. Contains unsigned "
                      "32-bit integers corresponding to mutations in the "
                      "population. (read-only)")
        .def("as_dict",
             // This lambda shows that
             // we can return dicts
             // with a mix of scalars + containers
             [](const fwdpp::gamete &g) noexcept {
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
        .def(py::pickle(
            [](const fwdpp::gamete &g) {
                return py::make_tuple(g.n, g.mutations, g.smutations);
            },
            [](py::tuple t) {
                return fwdpp::gamete(t[0].cast<fwdpp::uint_t>(),
                                     t[1].cast<std::vector<fwdpp::uint_t>>(),
                                     t[2].cast<std::vector<fwdpp::uint_t>>());
            }))
        .def("__eq__", [](const fwdpp::gamete &a, const fwdpp::gamete &b) {
            return a == b;
        });

    // Sugar types
    py::class_<fwdpp::popgenmut, fwdpp::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init<double, double, double, unsigned, std::uint16_t>(),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"),
             py::arg("label"),
             R"delim(
                Construct a mutations.

                :param pos: Mutation position (float)
                :param s: Effect size (float)
                :param h: Dominance term (float)
                :param g: Origin time (unsigned integer)
                :param label: Label (16 bit integer)

                .. testcode::

                    import fwdpy11
                    m = fwdpy11.Mutation(1.0, -1.0, 0.25, 0, 0)
                    print(m.pos)
                    print(m.s)
                    print(m.h)
                    print(m.g)
                    print(m.label)

                .. testoutput::
                    
                    1.0
                    -1.0
                    0.25
                    0
                    0
                )delim")
        .def(py::init<fwdpp::popgenmut::constructor_tuple>(),
             R"delim(
                Construct mutation from a tuple.

                The tuple should contain (pos, s, h, g, label)

                .. testcode::

                    import fwdpy11
                    m = fwdpy11.Mutation((1.0, -1.0, 0.25, 0, 0))
                    print(m.pos)
                    print(m.s)
                    print(m.h)
                    print(m.g)
                    print(m.label)

                .. testoutput::
                    
                    1.0
                    -1.0
                    0.25
                    0
                    0
                )delim")
        .def_readonly(
            "g", &fwdpp::popgenmut::g,
            "Generation when mutation arose (origination time). (read-only)")
        .def_readonly("s", &fwdpp::popgenmut::s,
                      "Selection coefficient/effect size. (read-only)")
        .def_readonly("h", &fwdpp::popgenmut::h,
                      "Dominance/effect in heterozygotes. (read-only)")
        .def_property_readonly(
            "key",
            [](const fwdpp::popgenmut &m) {
                return py::make_tuple(m.pos, m.s, m.g);
            },
            R"delim(It is often useful to have a unique key for
                    tracking mutations.  This property returns 
                    the tuple (pos, esize, origin).

                    .. versionadded:: 0.1.3.a1
                   )delim")
        .def(py::pickle(
            [](const fwdpp::popgenmut &m) {
                return py::make_tuple(m.pos, m.s, m.h, m.g, m.xtra);
            },
            [](py::tuple p) {
                return std::unique_ptr<fwdpp::popgenmut>(new fwdpp::popgenmut(
                    p[0].cast<double>(), p[1].cast<double>(),
                    p[2].cast<double>(), p[3].cast<unsigned>(),
                    p[4].cast<std::uint16_t>()));
            }))
        .def("__str__",
             [](const fwdpp::popgenmut &m) {
                 return "Mutation[" + std::to_string(m.pos) + ","
                        + std::to_string(m.s) + "," + std::to_string(m.h) + ","
                        + std::to_string(m.g) + "," + std::to_string(m.xtra)
                        + "]";
             })
        .def("__eq__", [](const fwdpp::popgenmut &a,
                          const fwdpp::popgenmut &b) { return a == b; });
}
