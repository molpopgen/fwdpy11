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
#include <pybind11/numpy.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpy11/types.hpp>

namespace py = pybind11;

static_assert(sizeof(char) == sizeof(std::int8_t),
              "sizeof(char) must equal sizeof(std::int8_t)");

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
            :param removeFixed: (boolean, defaults to True) Whether or not to include fixations.

            :rtype: tuple

            :return: A tuple.  The first element contains neutral variants, and the second
            contains selected variants.

            )delim",
          py::arg("rng"), py::arg("pop"), py::arg("samplesize"),
          py::arg("removeFixed") = true);

    m.def("sample_separate",
          [](const fwdpy11::singlepop_t &pop,
             const std::vector<std::size_t> &individuals,
             const bool removeFixed) {
              if (std::any_of(
                      individuals.begin(), individuals.end(),
                      [&pop](const std::size_t i) { return i >= pop.N; }))
                  {
                      throw std::runtime_error("individual index >= pop.N");
                  }
              return KTfwd::sample_separate(pop, individuals, removeFixed);
          },
          R"delim(
			Get a sample from a populaton based on a specific set of
            diploids.

            :param pop: A `fwdpy11.fwdpy11_types.Spop`
            :param individuals: (list of integers) The individuals to include in the sample. 
            :param removeFixed: (boolean, defaults to True) Whether or not to include fixations.

            :rtype: tuple

            :return: A tuple.  The first element contains neutral variants, and the second
            contains selected variants.

            )delim",
          py::arg("pop"), py::arg("individuals"),
          py::arg("removeFixed") = true);

    m.def("sample_separate",
          [](const fwdpy11::GSLrng_t &rng, const fwdpy11::multilocus_t &pop,
             const unsigned nsam, const bool removeFixed,
             const std::vector<std::pair<double, double>> &locus_boundaries) {
              return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed,
                                            locus_boundaries);
          });

    m.def("sample_separate",
          [](const fwdpy11::multilocus_t &pop,
             const std::vector<std::size_t> &individuals,
             const bool removeFixed,
             const std::vector<std::pair<double, double>> &locus_boundaries) {
              return KTfwd::sample_separate(pop, individuals, removeFixed,
                                            locus_boundaries);
          });

    py::class_<KTfwd::data_matrix>(m, "DataMatrix",
                                   R"delim(
		Represent a sample from a population in a matrix format.

		There are two possible representations of the data:

		1. As a genotype matrix, where individuals are encoded a 0,1, or 2
		copies of the derived mutation. There is one row per diploid here.
		2. As a haplotype matrix, with two rows per diploid, and each
		column containing a 0 (ancestral) or 1 (derived) label.

		You do not create objects of this type directly.  Instead, you use one of 
		the following functions:

		* :func:`fwdpy11.sampling.genotype_matrix`
		* :func:`fwdpy11.sampling.haplotype_matrix`
		
		Please see the examples in the fwdpy11 manual for examples of generating
		instances of this type.  The API requires multiple steps, in order to 
		maximize flexibility.
		)delim")
        .def(py::init<>())
        .def(py::init<std::size_t>())
        .def("neutral",
             [](const KTfwd::data_matrix &m) {
                 return py::array_t<std::int8_t>(
                     m.neutral.size(),
                     reinterpret_cast<const std::int8_t *>(m.neutral.data()));
             },
             "Return the neutral variants as a 1d numpy array.")
        .def("selected",
             [](const KTfwd::data_matrix &m) {
                 return py::array_t<std::int8_t>(
                     m.selected.size(),
                     reinterpret_cast<const std::int8_t *>(m.neutral.data()));
             },
             "Return the selected variants as a 1d numpy array.")
        .def_readonly("neutral_positions",
                      &KTfwd::data_matrix::neutral_positions,
                      "The list of neutral mutation positions.")
        .def_readonly("selected_positions",
                      &KTfwd::data_matrix::selected_positions,
                      "The list of selected mutation positions.")
        .def_readonly(
            "neutral_popfreq", &KTfwd::data_matrix::neutral_popfreq,
            "The list of population frequencies of neutral mutations.")
        .def_readonly(
            "selected_popfreq", &KTfwd::data_matrix::selected_popfreq,
            "The list of population frequencies of selected mutations.")
        .def_readonly("nrow", &KTfwd::data_matrix::nrow,
                      "Number of rows in the matrix.")
        .def("ncol_neutral",
             [](const KTfwd::data_matrix &m) {
                 return m.neutral.size() / m.nrow;
             },
             "Number of columns in the neutral variant matrix.")
        .def("ncol_selected",
             [](const KTfwd::data_matrix &m) {
                 return m.selected.size() / m.nrow;
             },
             "Number of columns in the selected variant matrix.")
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

#define MUTATION_KEYS(POPTYPE, CLASSTYPE)                                     \
    m.def("mutation_keys",                                                    \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const bool neutral, const bool selected) {                       \
              return KTfwd::mutation_keys(pop, individuals, neutral,          \
                                          selected);                          \
          },                                                                  \
          "Generate a tuple of (mutation_key, sample count) for mutations\n"  \
          "in a specific set of diploids in a :class:`" CLASSTYPE "` "        \
          "object.\n\n"                                                       \
          ":param pop: A population object.\n"                                \
          ":param individuals: A list of indexes to diploids.\n"              \
          ":param netural: (True) A boolean indicating whether to include "   \
          "neutral variants.\n"                                               \
          ":param selected: (True) A boolean indicating whether to include "  \
          "selected variants.\n\n"                                            \
          ".. note:: Mutation keys are unsorted.\n",                          \
          py::arg("pop"), py::arg("individuals"), py::arg("neutral") = true,  \
          py::arg("selected") = true);

    using keytype = std::vector<std::pair<std::size_t, KTfwd::uint_t>>;

#define GENOTYPE_MATRIX(POPTYPE, CLASSTYPE)                                   \
    m.def("genotype_matrix",                                                  \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const keytype &neutral_keys, const keytype &selected_keys) {     \
              return KTfwd::genotype_matrix<POPTYPE>(                         \
                  pop, individuals, neutral_keys, selected_keys);             \
          },                                                                  \
          "Generate a :class:fwdpy11.sampling.DataMatrix from a "             \
          ":class:`" CLASSTYPE "` object.\n"                                  \
          "The DataMatrix will be encoded as diploid genotypes.\n\n"          \
          ":param pop: A population object.\n"                                \
          ":param keys: The return value from "                               \
          ":func:`fwdpy11.sampling.mutation_keys`.\n\n"                       \
          ":rtype: :class:`fwdpy11.sampling.DataMatrix` encoded as a "        \
          "genotype matrix\n");

#define HAPLOTYPE_MATRIX(POPTYPE, CLASSTYPE)                                  \
    m.def("haplotype_matrix",                                                 \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const keytype &neutral_keys, const keytype &selected_keys) {     \
              return KTfwd::haplotype_matrix<POPTYPE>(                        \
                  pop, individuals, neutral_keys, selected_keys);             \
          },                                                                  \
          "Generate a :class:fwdpy11.sampling.DataMatrix from a "             \
          ":class:`" CLASSTYPE "` object.\n"                                  \
          "The DataMatrix will be encoded as haplotypes.\n\n"                 \
          ":param pop: A population object.\n"                                \
          ":param neutral_keys: The return value from "                       \
          ":func:`fwdpy11.sampling.mutation_keys`.\n\n"                       \
          ":param selected_keys: The return value from "                      \
          ":func:`fwdpy11.sampling.mutation_keys`.\n\n"                       \
          ":rtype: :class:`fwdpy11.sampling.DataMatrix` encoded as a "        \
          "haplotype matrix\n");

    MUTATION_KEYS(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.Spop");
    MUTATION_KEYS(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    MUTATION_KEYS(fwdpy11::singlepop_gm_vec_t,
                  "fwdpy11.fwdpy11_types.SpopGeneralMutVec");

    GENOTYPE_MATRIX(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.Spop");
    GENOTYPE_MATRIX(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    GENOTYPE_MATRIX(fwdpy11::singlepop_gm_vec_t,
                    "fwdpy11.fwdpy11_types.SpopGeneralMutVec");

    HAPLOTYPE_MATRIX(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.Spop");
    HAPLOTYPE_MATRIX(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    HAPLOTYPE_MATRIX(fwdpy11::singlepop_gm_vec_t,
                     "fwdpy11.fwdpy11_types.SpopGeneralMutVec");
    return m.ptr();
}
