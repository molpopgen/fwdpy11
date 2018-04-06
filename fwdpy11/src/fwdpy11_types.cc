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
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/fwdpp/sugar/sampling.hpp>

namespace py = pybind11;

PYBIND11_MODULE(fwdpy11_types, m)
{
    m.doc() = "Wrap C++ types specific to fwdpy11.";

    py::class_<fwdpy11::Diploid>(
        m, "SingleLocusDiploid",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_static("create", &fwdpy11::Diploid::create)
        .def_readonly("first", &fwdpy11::Diploid::first,
                      "Key to first gamete. (read-only)")
        .def_readonly("second", &fwdpy11::Diploid::second,
                      "Key to second gamete. (read-only)")
        .def_readonly("w", &fwdpy11::Diploid::w, "Fitness. (read-only)")
        .def_readonly("g", &fwdpy11::Diploid::g, "Genetic value (read-only).")
        .def_readonly("e", &fwdpy11::Diploid::e,
                      "Random/environmental effects (read-only).")
        .def_readonly("label", &fwdpy11::Diploid::label,
                      "Index of the diploid in its deme")
        .def_readonly("deme", &fwdpy11::Diploid::deme,
                      R"delim(
                Deme label for individual.

                .. versionadded:: 0.1.5
                )delim")
        .def_readonly("sex", &fwdpy11::Diploid::sex,
                      R"delim(
                Sex label for individual.

                .. versionadded:: 0.1.5
                )delim")
        .def_readonly("parental_data", &fwdpy11::Diploid::parental_data,
                      R"delim(
				Python object representing information about parents.
				The details are simulation-dependent.

				.. versionadded:: 0.1.4
				)delim")
        .def(py::pickle(
            [](const fwdpy11::Diploid& d) {
                return py::make_tuple(d.first, d.second, d.w, d.g, d.e,
                                      d.label, d.parental_data, d.deme, d.sex);
            },
            [](py::tuple t) {
                std::unique_ptr<fwdpy11::Diploid> d(new fwdpy11::Diploid(
                    t[0].cast<std::size_t>(), t[1].cast<std::size_t>()));
                d->w = t[2].cast<double>();
                d->g = t[3].cast<double>();
                d->e = t[4].cast<double>();
                d->label = t[5].cast<decltype(fwdpy11::Diploid::label)>();
                d->parental_data
                    = t[6].cast<std::tuple<std::size_t, std::size_t>>();
                d->deme = t[7].cast<std::uint32_t>();
                d->sex = t[8].cast<std::int32_t>();
                return d;
            }))
        .def("__eq__", [](const fwdpy11::Diploid& a,
                          const fwdpy11::Diploid& b) { return a == b; });
}
