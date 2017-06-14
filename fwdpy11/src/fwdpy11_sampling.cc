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

/* This file exposes the same fwdpp functions over and over again
 * for different population types.  We use macros to reduce redundancy.
 */

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <array>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/internal/IOhelp.hpp>
#include <fwdpy11/types.hpp>
#include <gsl/gsl_matrix_char.h>
namespace py = pybind11;

static_assert(sizeof(char) == sizeof(std::int8_t),
              "sizeof(char) must equal sizeof(std::int8_t)");

py::list
matrix_to_sample(const std::vector<std::int8_t> &data,
                 const std::vector<double> &pos, const std::size_t nrow)
// returns a data structure compatible with libsequence/pylibseq iff
// the data correspond to a haplotype matrix
{
    std::size_t ncol = data.size() / nrow;
    const std::array<std::int8_t, 3> states{ '0', '1', '2' };
    auto v = gsl_matrix_char_const_view_array(
        reinterpret_cast<const char *>(data.data()), nrow, ncol);
    py::list rv;
    for (std::size_t i = 0; i < ncol; ++i)
        {
            auto c = gsl_matrix_char_const_column(&v.matrix, i);
            std::string column_data;
            for (std::size_t j = 0; j < c.vector.size; ++j)
                {
                    column_data.push_back(states[static_cast<std::int8_t>(
                        gsl_vector_char_get(&c.vector, j))]);
                }
            if (column_data.size() != nrow)
                {
                    throw std::runtime_error("column_data.size() != nrow");
                }
            rv.append(py::make_tuple(pos[i], std::move(column_data)));
        }
    return rv;
}

py::dict
separate_samples_by_loci(
    const std::vector<std::pair<double, double>> &boundaries, py::list sample)
// For a multi-locus pop, it is convenient to split samples by
// loci.  This function does that using pop.locus_boundaries.
// If pop.locus_boundaries is not properly set, an exception
// is likely going to be triggered
{
    py::dict rv;
    if (sample.size() == 0)
        {
            return rv;
        }
    for (std::size_t i = 0; i < boundaries.size(); ++i)
        {
            rv[py::int_(i)] = py::list();
        }
    for (auto &&item : sample)
        {
            py::tuple site = py::reinterpret_borrow<py::tuple>(item);
            if (site.size() != 2)
                {
                    throw std::runtime_error("invalid tuple length: "
                                             + std::to_string(site.size())
                                             + " seen when 2 was expected");
                }
            auto itr
                = std::find_if(boundaries.begin(), boundaries.end(),
                               [&site](const std::pair<double, double> &b) {
                                   return site[0].cast<double>() >= b.first
                                          && site[0].cast<double>() < b.second;
                               });
            if (itr == boundaries.end())
                {
                    throw std::runtime_error(
                        "could not find locus for mutation at position"
                        + std::to_string(site[0].cast<double>()));
                }
            auto d = std::distance(boundaries.begin(), itr);
            py::list li = py::reinterpret_borrow<py::list>(rv[py::int_(d)]);
            li.append(site);
        }
    return rv;
}

PYBIND11_MAKE_OPAQUE(std::vector<std::int8_t>);

PYBIND11_PLUGIN(sampling)
{
    py::module m("sampling", "Taking samples from populations");

#define SAMPLE_SEPARATE_RANDOM(POPTYPE, CLASSTYPE)                            \
    m.def("sample_separate",                                                  \
          [](const fwdpy11::GSLrng_t &rng, const POPTYPE &pop,                \
             const KTfwd::uint_t samplesize, const bool removeFixed) {        \
              return KTfwd::sample_separate(rng.get(), pop, samplesize,       \
                                            removeFixed);                     \
          },                                                                  \
          "Take a sample of :math:`n` chromosomes from a population "         \
          "(`n/2` diploids.\n\n"                                              \
          ":param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`\n"             \
          ":param pop: A :class:`" CLASSTYPE "`\n"                            \
          ":param samplesize: (int) The sample size.\n"                       \
          ":param removeFixed: (boolean, defaults to True) Whether or not to" \
          "include fixations.\n"                                              \
          ":rtype: tuple\n\n"                                                 \
          ":return: A tuple.  The first element contains neutral variants,"   \
          "and the second contains selected variants.\n",                     \
          py::arg("rng"), py::arg("pop"), py::arg("samplesize"),              \
          py::arg("removeFixed") = true);

#define SAMPLE_SEPARATE_IND(POPTYPE, CLASSTYPE)                               \
    m.def("sample_separate",                                                  \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const bool removeFixed) {                                        \
              return KTfwd::sample_separate(pop, individuals, removeFixed);   \
          },                                                                  \
          "Take a sample of specific individuals from a population.\n\n "     \
          ":param pop: A :class:`" CLASSTYPE "`\n"                            \
          ":param individuals : (list of int) Indexes of individuals.\n"      \
          ":param removeFixed: (boolean, defaults to True) Whether or not to" \
          "include fixations.\n"                                              \
          ":rtype: tuple\n\n"                                                 \
          ":return: A tuple.  The first element contains neutral variants,"   \
          "and the second contains selected variants.\n",                     \
          py::arg("pop"), py::arg("individuals"),                             \
          py::arg("removeFixed") = true);

    SAMPLE_SEPARATE_RANDOM(fwdpy11::singlepop_t,
                           "fwdpy11.fwdpy11_types.SlocusPop")
    SAMPLE_SEPARATE_RANDOM(fwdpy11::multilocus_t,
                           "fwdpy11.fwdpy11_types.MlocusPop")
    SAMPLE_SEPARATE_RANDOM(fwdpy11::singlepop_gm_vec_t,
                           "fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec")
    SAMPLE_SEPARATE_IND(fwdpy11::singlepop_t,
                        "fwdpy11.fwdpy11_types.SlocusPop")
    SAMPLE_SEPARATE_IND(fwdpy11::multilocus_t,
                        "fwdpy11.fwdpy11_types.MlocusPop")
    SAMPLE_SEPARATE_IND(fwdpy11::singlepop_gm_vec_t,
                        "fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec")

    py::bind_vector<std::vector<std::int8_t>>(m, "Vec8",
                                              py::buffer_protocol());

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
		
		Please see :ref:`datamatrix` for examples of generating
		instances of this type.  The API requires multiple steps, in order to 
		maximize flexibility.
		)delim")
        .def(py::init<>())
        .def(py::init<std::size_t>())
        .def_readonly("neutral", &KTfwd::data_matrix::neutral,
                      R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array
                )delim")
        .def_readonly("selected", &KTfwd::data_matrix::selected,
                      R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array
                )delim")
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
        .def("ndim_neutral",
             [](const KTfwd::data_matrix &dm) {
                 return py::make_tuple(dm.nrow, dm.neutral.size() / dm.nrow);
             },
             R"delim(
             Return the dimensions of the neutral matrix
             
             :rtype: tuple

             .. versionadded:: 0.1.2
                Replaces ncol and nrow_neutral functions
             )delim")
        .def("ndim_selected",
             [](const KTfwd::data_matrix &dm) {
                 return py::make_tuple(dm.nrow, dm.selected.size() / dm.nrow);
             },
             R"delim(
             Return the dimensions of the selected matrix

             :rtype: tuple
             
             .. versionadded:: 0.1.2
                Replaces ncol and nrow_selected functions
             )delim")
        .def("__getstate__",
             [](const KTfwd::data_matrix &d) {
                 std::ostringstream o;
                 KTfwd::fwdpp_internal::scalar_writer w;
                 w(o, &d.nrow, 1);
                 auto nsites = d.neutral_positions.size();
                 w(o, &nsites, 1);
                 if (nsites)
                     {
                         auto l = d.neutral.size();
                         w(o, &l);
                         w(o, d.neutral.data(), d.neutral.size());
                         w(o, d.neutral_positions.data(), nsites);
                         w(o, d.neutral_popfreq.data(), nsites);
                     }
                 nsites = d.selected_positions.size();
                 w(o, &nsites, 1);
                 if (nsites)
                     {
                         auto l = d.neutral.size();
                         w(o, &l);
                         w(o, d.selected.data(), d.selected.size());
                         w(o, d.selected_positions.data(), nsites);
                         w(o, d.selected_popfreq.data(), nsites);
                     }
                 return py::bytes(o.str());
             })
        .def("__setstate__", [](KTfwd::data_matrix &d, py::bytes b) {
            std::istringstream data(b);
            KTfwd::fwdpp_internal::scalar_reader r;
            std::size_t n, n2;
            r(data, &n);
            new (&d) KTfwd::data_matrix(n);
            r(data, &n);
            if (n)
                {
                    r(data, &n2);
                    d.neutral.resize(n2);
                    r(data, d.neutral.data(), n2);
                    d.neutral_positions.resize(n);
                    r(data, d.neutral_positions.data(), n);
                    d.neutral_popfreq.resize(n);
                    r(data, d.neutral_popfreq.data(), n);
                }
            r(data, &n);
            if (n)
                {
                    r(data, &n2);
                    d.selected.resize(n2);
                    r(data, d.selected.data(), n2);
                    d.selected_positions.resize(n);
                    r(data, d.selected_positions.data(), n);
                    d.selected_popfreq.resize(n);
                    r(data, d.selected_popfreq.data(), n);
                }
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

    MUTATION_KEYS(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.SlocusPop");
    MUTATION_KEYS(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    MUTATION_KEYS(fwdpy11::singlepop_gm_vec_t,
                  "fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec");

    GENOTYPE_MATRIX(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.SlocusPop");
    GENOTYPE_MATRIX(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    GENOTYPE_MATRIX(fwdpy11::singlepop_gm_vec_t,
                    "fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec");

    HAPLOTYPE_MATRIX(fwdpy11::singlepop_t, "fwdpy11.fwdpy11_types.SlocusPop");
    HAPLOTYPE_MATRIX(fwdpy11::multilocus_t, "fwdpy11.fwdpy11_types.MlocusPop");
    HAPLOTYPE_MATRIX(fwdpy11::singlepop_gm_vec_t,
                     "fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec");

    m.def("matrix_to_sample",
          [](const KTfwd::data_matrix &m, const bool neutral)

          {
              return (neutral) ? matrix_to_sample(m.neutral,
                                                  m.neutral_positions, m.nrow)
                               : matrix_to_sample(
                                     m.selected, m.selected_positions, m.nrow);
          },
          R"delim(
          Convert a :class:`fwdpy11.sampling.DataMatrix` into a
          list of tuples (float,string).

          .. versionadded:: 0.1.1

          :param m: A :class:`fwdpy11.sampling.DataMatrix`
          :param neutral: (True) Return data for neutral or selected sites.

          :rtype: list of tuples

          :return: The data in list format.  Positions are in same order 
              as the original matrix.  The genotypes are represented as strings
              with elements 0 through nrow-1 being in the same order as the original matrix.

          .. note::
            If the input data represent a haplotype matrix, then the return value
            may be further processed using `pylibseq <http://molpopgen.github.io/pylibseq/>`_.
            Any filtering on position, frequency, etc., should have aleady been done when the 
            original matrix was generated.  However, you may filter the return value as you see fit.
            If the matrix was created from a multi-locus simulation, you may wish to use
            :func:`fwdpy11.sampling.separate_samples_by_loci` to split the return value up
			into separate lists for each locus.

          )delim",
          py::arg("m"), py::arg("neutral") = true);

    m.def("separate_samples_by_loci", &separate_samples_by_loci,
          R"delim(
            Convert the output from :func:`fwdpy11.sampling.matrix_to_sample` into 
            separate records per locus.

			.. versionadded:: 0.1.1

			.. versionchanged:: 0.1.2
				Take a list of positions as arguments and not a population object.
		
            :param boundaries: A list of [start,stop) tuples representing positions.
            :param sample: The return value of :func:`fwdpy11.sampling.matrix_to_sample`

            :rtype: dict of lists of tuples

            :return: The data returned follow the same structure 
				as :func:`fwdpy11.sampling.matrix_to_sample`,
                but there is one entry per locus.  
				The key for each entry in the dict is the locus index.
            )delim");

    return m.ptr();
}
