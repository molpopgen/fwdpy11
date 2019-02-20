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
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpy11/types/Mutation.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_MODULE(fwdpy11_types, m)
{
    m.doc() = "Wrap C++ types specific to fwdpy11.";

    py::bind_vector<std::vector<double>>(
        m, "VecDouble", "C++ vector of 64 bit floats.", py::buffer_protocol())
        .def(py::pickle(
            [](const std::vector<double> &v) {
                py::list rv;
                for (auto &&i : v)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<double> rv;
                for (auto &&i : l)
                    {
                        rv.push_back(i.cast<double>());
                    }
                return rv;
            }));

    // Sugar types
    py::class_<fwdpy11::Mutation, fwdpp::mutation_base>(
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
        .def(py::init([](double pos, double s, double h, fwdpp::uint_t g,
                         py::list esizes, py::list heffects,
                         std::uint16_t label) {
                 std::vector<double> esizes_;
                 std::vector<double> heffects_;
                 for (auto i : esizes)
                     {
                         esizes_.push_back(i.cast<double>());
                     }
                 for (auto i : heffects)
                     {
                         heffects_.push_back(i.cast<double>());
                     }
                 return fwdpy11::Mutation(pos, s, h, g, std::move(esizes_),
                                          std::move(heffects_), label);
             }),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"),
             py::arg("esizes"), py::arg("heffects"), py::arg("label"),
             R"delim(
			 Construct a mutation with both a constant effect and/or
			 a vector of effects.

			 :param pos: Mutation position (float)
			 :param s: Effect size (float)
			 :param h: Dominance term (float)
			 :param g: Origin time (unsigned integer)
			 :param esizes: List of effect sizes (list of float)
			 :param heffects: List of heterozygous effects (list of float)
			 :param label: Label (16 bit integer)
				
			 .. versionadded:: 0.2.0

			 )delim")
        .def(py::init<fwdpy11::Mutation::constructor_tuple>(),
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
            "g", &fwdpy11::Mutation::g,
            "Generation when mutation arose (origination time). (read-only)")
        .def_readonly("s", &fwdpy11::Mutation::s,
                      "Selection coefficient/effect size. (read-only)")
        .def_readonly("h", &fwdpy11::Mutation::h,
                      "Dominance/effect in heterozygotes. (read-only)")
        .def_readonly("heffects", &fwdpy11::Mutation::heffects,
                      R"delim(
				Vector of heterozygous effects.

				.. versionadded:: 0.2.0
				)delim")
        .def_readonly("esizes", &fwdpy11::Mutation::esizes,
                      R"delim(
				Vector of effect sizes.

				.. versionadded:: 0.2.0
				)delim")
        .def_property_readonly(
            "key",
            [](const fwdpy11::Mutation &m) {
                return py::make_tuple(m.pos, m.s, m.g);
            },
            R"delim(It is often useful to have a unique key for
                    tracking mutations.  This property returns 
                    the tuple (pos, esize, origin).

                    .. versionadded:: 0.1.3.a1
                   )delim")
        .def(py::pickle(
            [](const fwdpy11::Mutation &m) {
                return py::make_tuple(m.pos, m.s, m.h, m.g, m.esizes,
                                      m.heffects, m.xtra);
            },
            [](py::tuple p) {
                return std::unique_ptr<fwdpy11::Mutation>(
                    new fwdpy11::Mutation(
                        p[0].cast<double>(), p[1].cast<double>(),
                        p[2].cast<double>(), p[3].cast<unsigned>(),
                        p[4].cast<std::vector<double>>(),
                        p[5].cast<std::vector<double>>(),
                        p[6].cast<std::uint16_t>()));
            }))
        .def("__str__",
             [](const fwdpy11::Mutation &m) {
                 return "Mutation[" + std::to_string(m.pos) + ","
                        + std::to_string(m.s) + "," + std::to_string(m.h) + ","
                        + std::to_string(m.g) + "," + std::to_string(m.xtra)
                        + "]";
             })
        .def("__eq__", [](const fwdpy11::Mutation &a,
                          const fwdpy11::Mutation &b) { return a == b; });

    py::class_<fwdpy11::DiploidGenotype>(
        m, "DiploidGenotype",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_readonly("first", &fwdpy11::DiploidGenotype::first,
                      "Key to first gamete. (read-only)")
        .def_readonly("second", &fwdpy11::DiploidGenotype::second,
                      "Key to second gamete. (read-only)")
        .def(py::pickle(
            [](const fwdpy11::DiploidGenotype &d) {
                return py::make_tuple(d.first, d.second);
            },
            [](py::tuple t) {
                std::unique_ptr<fwdpy11::DiploidGenotype> d(
                    new fwdpy11::DiploidGenotype{ t[0].cast<std::size_t>(),
                                                  t[1].cast<std::size_t>() });
                return d;
            }))
        .def("__eq__",
             [](const fwdpy11::DiploidGenotype &a,
                const fwdpy11::DiploidGenotype &b) { return a == b; });

    // TODO: pickling support
    py::class_<fwdpy11::DiploidMetadata>(m, "DiploidMetadata",
                                         "Diploid meta data.")
        .def_readwrite("g", &fwdpy11::DiploidMetadata::g, "Genetic value.")
        .def_readwrite("e", &fwdpy11::DiploidMetadata::e,
                       "Random component of trait value.")
        .def_readwrite("w", &fwdpy11::DiploidMetadata::w, "Fitness.")
        .def_property(
            "geography",
            [](const fwdpy11::DiploidMetadata &d) {
                return py::make_tuple(d.geography[0], d.geography[1],
                                      d.geography[2]);
            },
            [](fwdpy11::DiploidMetadata &d,
               const std::tuple<double, double, double> &input) {
                d.geography[0] = std::get<0>(input);
                d.geography[1] = std::get<1>(input);
                d.geography[2] = std::get<2>(input);
            },
            "Array containing the geographic location of the individual.")
        .def_property("parents",
                      [](const fwdpy11::DiploidMetadata &d) {
                          return py::make_tuple(d.parents[0], d.parents[1]);
                      },
                      [](fwdpy11::DiploidMetadata &d,
                         const std::pair<std::size_t, std::size_t> &input) {
                          d.parents[0] = input.first;
                          d.parents[1] = input.second;
                      },
                      "Array containing the label fields of the parents.")
        .def_readwrite("sex", &fwdpy11::DiploidMetadata::sex, "Sex.")
        .def_readwrite("deme", &fwdpy11::DiploidMetadata::deme, "Deme.")
        .def_readwrite("label", &fwdpy11::DiploidMetadata::label,
                       "Index of the individual in the population.");
}
