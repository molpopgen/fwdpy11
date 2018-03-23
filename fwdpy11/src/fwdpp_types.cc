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
//#include <fwdpp/sugar/generalmut.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::uint_t>);
//PYBIND11_MAKE_OPAQUE(
//    std::vector<double>); // for generalmut_vec::s and generalmut_vec::h

PYBIND11_MODULE(fwdpp_types, m)
{
    m.doc() = "Wrap C++ types from fwdpp.";

    // low-level types

    py::class_<KTfwd::mutation_base>(m, "MutationBase",
                                     R"delim(
                                        Base class for mutations.
                                     )delim")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readonly("pos", &KTfwd::mutation_base::pos, "Position (float).")
        .def_readonly("neutral", &KTfwd::mutation_base::neutral, "Boolean")
        .def_readwrite("label", &KTfwd::mutation_base::xtra,
                       "A 16-bit unsigned integer that can be used for adding "
                       "\"meta-data\" to mutations");

    py::bind_vector<std::vector<KTfwd::uint_t>>(
        m, "VecUint32", "Vector of unsigned 32-bit integers.",
        py::buffer_protocol())
        .def(py::pickle(
            [](const std::vector<KTfwd::uint_t>& v) {
                py::list rv;
                for (auto&& i : v)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<KTfwd::uint_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<KTfwd::uint_t>());
                    }
                return rv;
            }));

    py::class_<KTfwd::gamete>(m, "Gamete", R"delim(
    A gamete.  This object represents a haplotype
    in a contiguous genomic region.
)delim")
        .def(py::init<KTfwd::gamete::constructor_tuple>(),
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
        .def_readonly("n", &KTfwd::gamete::n,
                      "Number of occurrences in the population. This has "
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
        .def(py::pickle(
            [](const KTfwd::gamete &g) {
                return py::make_tuple(g.n, g.mutations, g.smutations);
            },
            [](py::tuple t) {
                return KTfwd::gamete(t[0].cast<KTfwd::uint_t>(),
                                     t[1].cast<std::vector<KTfwd::uint_t>>(),
                                     t[2].cast<std::vector<KTfwd::uint_t>>());
            }))
        .def("__eq__", [](const KTfwd::gamete &a, const KTfwd::gamete &b) {
            return a == b;
        });

    // Sugar types
    py::class_<KTfwd::popgenmut, KTfwd::mutation_base>(
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
        .def(py::init<KTfwd::popgenmut::constructor_tuple>(),
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
            "g", &KTfwd::popgenmut::g,
            "Generation when mutation arose (origination time). (read-only)")
        .def_readonly("s", &KTfwd::popgenmut::s,
                      "Selection coefficient/effect size. (read-only)")
        .def_readonly("h", &KTfwd::popgenmut::h,
                      "Dominance/effect in heterozygotes. (read-only)")
        .def_property_readonly(
            "key",
            [](const KTfwd::popgenmut &m) {
                return py::make_tuple(m.pos, m.s, m.g);
            },
            R"delim(It is often useful to have a unique key for
                    tracking mutations.  This property returns 
                    the tuple (pos, esize, origin).

                    .. versionadded:: 0.1.3.a1
                   )delim")
        .def(py::pickle(
            [](const KTfwd::popgenmut &m) {
                return py::make_tuple(m.pos, m.s, m.h, m.g, m.xtra);
            },
            [](py::tuple p) {
                return std::unique_ptr<KTfwd::popgenmut>(new KTfwd::popgenmut(
                    p[0].cast<double>(), p[1].cast<double>(),
                    p[2].cast<double>(), p[3].cast<unsigned>(),
                    p[4].cast<std::uint16_t>()));
            }))
        .def("__str__",
             [](const KTfwd::popgenmut &m) {
                 return "Mutation[" + std::to_string(m.pos) + ","
                        + std::to_string(m.s) + "," + std::to_string(m.h) + ","
                        + std::to_string(m.g) + "," + std::to_string(m.xtra)
                        + "]";
             })
        .def("__eq__", [](const KTfwd::popgenmut &a,
                          const KTfwd::popgenmut &b) { return a == b; });

    //py::bind_vector<std::vector<double>>(
    //    m, "VecDouble", "Vector of 64-bit floats.", py::buffer_protocol());

    //py::class_<KTfwd::generalmut_vec, KTfwd::mutation_base>(
    //    m, "GeneralMutVec",
    //    "Mutation type with vector of effect size and dominance terms.")
    //    .def(py::init<KTfwd::generalmut_vec::constructor_tuple>(),
    //         R"delim(
    //            Construct from a tuple.
    //            
    //            .. testcode::
    //                
    //                import fwdpy11
    //                s = fwdpy11.VecDouble([-0.1, 0.1])
    //                h = fwdpy11.VecDouble([0.0, 1.0])
    //                m = fwdpy11.GeneralMutVec((s, h, 0.1, 3, 0))
    //                print(m.pos)
    //                print(m.g)
    //                print(m.label)
    //                print(list(m.s))
    //                print(list(m.h))

    //            .. testoutput::

    //                0.1
    //                3
    //                0
    //                [-0.1, 0.1]
    //                [0.0, 1.0]
    //            )delim")
    //    .def_readonly("s", &KTfwd::generalmut_vec::s,
    //                  "List of selection coefficients/effect sizes.")
    //    .def_readonly("h", &KTfwd::generalmut_vec::h,
    //                  "List of dominance terms.")
    //    .def_readonly("g", &KTfwd::generalmut_vec::g,
    //                  "Generation when mutation arose.")
    //    .def(py::pickle(
    //        [](const KTfwd::generalmut_vec &m) {
    //            py::list h, s;
    //            for (auto &&i : m.h)
    //                {
    //                    h.append(i);
    //                }
    //            for (auto &&i : m.s)
    //                {
    //                    s.append(i);
    //                }
    //            return py::make_tuple(m.pos, s, h, m.g, m.xtra);
    //        },
    //        [](py::tuple t) {
    //            py::list s = t[1].cast<py::list>();
    //            py::list h = t[2].cast<py::list>();
    //            std::vector<double> vs, vh;
    //            for (auto &&i : s)
    //                {
    //                    vs.push_back(i.cast<double>());
    //                }
    //            for (auto &&i : h)
    //                {
    //                    vh.push_back(i.cast<double>());
    //                }
    //            double pos = t[0].cast<double>();
    //            KTfwd::uint_t g = t[3].cast<KTfwd::uint_t>();
    //            auto xtra = t[4].cast<decltype(KTfwd::generalmut_vec::xtra)>();
    //            return KTfwd::generalmut_vec(std::make_tuple(
    //                std::move(vs), std::move(vh), pos, g, xtra));
    //        }))
    //    .def("__eq__", [](const KTfwd::generalmut_vec &a,
    //                      const KTfwd::generalmut_vec &b) { return a == b; });
}
