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
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpy11/types.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(sampling)
{
    py::module m("sampling", "Taking samples from populations");

    m.def("sample_separate",
          [](const fwdpy11::GSLrng_t &rng, const fwdpy11::singlepop_t &pop,
             const KTfwd::uint_t samplesize, const bool removeFixed) {
              return KTfwd::sample_separate(rng.get(), pop, samplesize,
                                            removeFixed);
          },
          R"delim(
            Take a sample of :math:`n` chromosomes from a population (`n/2` 
            diploids.

            :param rng: A `fwdpy11.fwdpy11_types.GSLrng`
            :param pop: A `fwdpy11.fwdpy11_types.Spop`
            :param samplesize: (int) The sample size.
            :param removeFixed (boolean, defaults to True) Whether or not to include fixations.

            :rtype: tuple

            :return: A tuple.  The first element contains neutral variants, and the second
            contains selected variants.

            )delim",
          py::arg("rng"), py::arg("pop"), py::arg("samplesize"),
          py::arg("removeFixed") = true);
    m.def(
        "sample_separate",
        [](const fwdpy11::GSLrng_t &rng, const fwdpy11::multilocus_t &pop,
           const unsigned nsam, const bool removeFixed,
           const std::vector<std::pair<double, double>> &locus_boundaries) {
            return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed,
                                          locus_boundaries);
        });

    py::class_<KTfwd::data_matrix>(m, "DataMatrix")
        .def(py::init<>())
        .def(py::init<std::size_t>())
        .def_readonly("neutral", &KTfwd::data_matrix::neutral)
        .def_readonly("selected", &KTfwd::data_matrix::selected)
        .def_readonly("neutral_positions",
                      &KTfwd::data_matrix::neutral_positions)
        .def_readonly("selected_positions",
                      &KTfwd::data_matrix::selected_positions)
        .def_readonly("neutral_popfreq", &KTfwd::data_matrix::neutral_popfreq)
        .def_readonly("selected_popfreq",
                      &KTfwd::data_matrix::selected_popfreq)
        .def_readonly("nrow", &KTfwd::data_matrix::nrow)
        .def("__getstate__",
             [](const KTfwd::data_matrix &d) {
                 return py::make_tuple(d.nrow, d.neutral, d.selected,
                                       d.neutral_positions,
                                       d.selected_positions, d.neutral_popfreq,
                                       d.selected_popfreq);
             })
        .def("__setstate__", [](KTfwd::data_matrix &d, py::tuple p) {
            new (&d) KTfwd::data_matrix(p[0].cast<std::size_t>());
            d.nrow = p[0].cast<std::size_t>();
            d.neutral = p[1].cast<std::vector<char>>();
            d.selected = p[1].cast<std::vector<char>>();
            d.neutral_positions = p[1].cast<std::vector<double>>();
            d.selected_positions = p[1].cast<std::vector<double>>();
            d.neutral_popfreq = p[1].cast<std::vector<double>>();
            d.selected_popfreq = p[1].cast<std::vector<double>>();
        });

    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy11::singlepop_t>);
    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy11::multilocus_t>);
    m.def("genotype_matrix", &KTfwd::genotype_matrix<fwdpy11::singlepop_t>);
    m.def("genotype_matrix", &KTfwd::genotype_matrix<fwdpy11::multilocus_t>);
    return m.ptr();
}
