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
#include <pybind11/numpy.h>
#include <array>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/data_matrix.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/sampling/data_matrix_functions.hpp>
#include <gsl/gsl_matrix_char.h>
namespace py = pybind11;

static_assert(sizeof(char) == sizeof(std::int8_t),
              "sizeof(char) must equal sizeof(std::int8_t)");

PYBIND11_MAKE_OPAQUE(std::vector<std::int8_t>);

PYBIND11_MODULE(sampling, m)
{
    m.doc() = "Taking samples from populations";

    py::class_<fwdpp::state_matrix>(m, "StateMatrix", py::buffer_protocol(),
            R"delim(
            Simple matrix representation of variation data.

            These are not constructed directly.  Rather,
            they are generated when a 
            :class:`fwdpy11.sampling.DataMatrix` is generated.

            This object supports the buffer protocol.

            .. versionadded:: 0.2.0
            )delim")
        .def_property_readonly(
            "shape",
            [](const fwdpp::state_matrix &sm) {
                if (sm.positions.empty())
                    {
                        return py::make_tuple(0, 0);
                    }
                if (sm.data.empty())
                    {
                        throw std::runtime_error("StatMatrix data are empty");
                    }
                return py::make_tuple(sm.positions.size(),
                                      sm.data.size() / sm.positions.size());
            },"Shape of the matrix.")
        .def_readonly("positions", &fwdpp::state_matrix::positions,
                      "The mutation positions")
        .def_buffer([](const fwdpp::state_matrix &sm) -> py::buffer_info {
            using value_type = std::int8_t;
            auto nrow = sm.positions.size();
            auto ncol = (nrow > 0) ? sm.data.size() / nrow : 0;
            return py::buffer_info(
                const_cast<value_type *>(sm.data.data()), sizeof(value_type),
                py::format_descriptor<value_type>::format(), 2, { nrow, ncol },
                { sizeof(value_type) * ncol, sizeof(value_type) });
        });

    py::class_<fwdpp::data_matrix>(m, "DataMatrix",
                                   R"delim(
		Represent a sample from a population in a matrix format.

		There are two possible representations of the data:

		1. As a genotype matrix, where individuals are encoded a 0,1, or 2
		copies of the derived mutation. There is one column per diploid here,
        and one row per variable site.

		2. As a haplotype matrix, with two columns per diploid, and each
		column containing a 0 (ancestral) or 1 (derived) label. Each row
        represents a variable site.

		You do not create objects of this type directly.  Instead, you use one of 
		the following functions:

		* :func:`fwdpy11.sampling.genotype_matrix`
		* :func:`fwdpy11.sampling.haplotype_matrix`
		
		Please see :ref:`datamatrix` for examples of generating
		instances of this type.  The API requires multiple steps, in order to 
		maximize flexibility.

        .. versionchanged:: 0.2.0

            Changed layout to row = variable site. 
            Changed to match fwdpp 0.7.0 layour where the neutral
            and selected data are represented as a 
            :class:`fwdpy11.sampling.StateMatrix`
		)delim")
        .def_readwrite("neutral", &fwdpp::data_matrix::neutral,
                       R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array

                .. versionchanged:: 0.1.4
                    Allow read/write access instead of readonly

                .. versionchanged:: 0.2.0
                    Type is :class:`fwdpy11.sampling.StateMatrix`
                )delim")
        .def_readwrite("selected", &fwdpp::data_matrix::selected,
                       R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array

                .. versionchanged:: 0.1.4
                    Allow read/write access instead of readonly

                .. versionchanged:: 0.2.0
                    Type is :class:`fwdpy11.sampling.StateMatrix`
                )delim")
        .def_readonly("ncol", &fwdpp::data_matrix::ncol,
                      "Sample size of the matrix")
        .def_readonly("neutral_keys", &fwdpp::data_matrix::neutral_keys,
                      "Keys for neutral mutations used to generate matrix")
        .def_readonly("selected_keys", &fwdpp::data_matrix::selected_keys,
                      "Keys for selected mutations used to generate matrix")
        .def(py::pickle(
            [](const fwdpp::data_matrix &d) {
                std::ostringstream o;
                fwdpp::io::scalar_writer w;
                auto nsites = d.neutral.positions.size();
                auto dsize = d.neutral.data.size();
                w(o, &nsites, 1);
                w(o, &dsize, 1);
                if (nsites)
                    {
                        w(o, d.neutral.data.data(), d.neutral.data.size());
                        w(o, d.neutral.positions.data(), nsites);
                        w(o, d.neutral_keys.data(), nsites);
                    }
                nsites = d.selected.positions.size();
                dsize = d.selected.data.size();
                w(o, &nsites, 1);
                w(o, &dsize, 1);
                if (nsites)
                    {
                        w(o, d.selected.data.data(), d.selected.data.size());
                        w(o, d.selected.positions.data(), nsites);
                        w(o, d.selected_keys.data(), nsites);
                    }
                return py::bytes(o.str());
            },
            [](py::bytes b) {
                std::istringstream data(b);
                fwdpp::io::scalar_reader r;
                std::size_t nsites, dsize;
                r(data, &nsites);
                r(data, &dsize);
                fwdpp::data_matrix d(dsize);
                if (nsites)
                    {
                        d.neutral.data.resize(dsize);
                        r(data, d.neutral.data.data(), dsize);
                        d.neutral.positions.resize(nsites);
                        r(data, d.neutral.positions.data(), nsites);
                        d.neutral_keys.resize(nsites);
                        r(data, d.neutral_keys.data(), nsites);
                    }
                r(data, &nsites);
                r(data, &dsize);
                if (nsites)
                    {
                        d.selected.data.resize(dsize);
                        r(data, d.selected.data.data(), dsize);
                        d.selected.positions.resize(nsites);
                        r(data, d.selected.positions.data(), nsites);
                        d.selected_keys.resize(nsites);
                        r(data, d.selected_keys.data(), nsites);
                    }
                return std::unique_ptr<fwdpp::data_matrix>(
                    new fwdpp::data_matrix(std::move(d)));
            }));

#define MUTATION_KEYS(POPTYPE, CLASSTYPE)                                     \
    m.def("mutation_keys",                                                    \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const bool neutral, const bool selected) {                       \
              return fwdpp::mutation_keys(pop, individuals, neutral,          \
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
          ".. note:: Mutation keys are returned unsorted.\n",                 \
          py::arg("pop"), py::arg("individuals"), py::arg("neutral") = true,  \
          py::arg("selected") = true);

    using keytype = std::vector<std::pair<std::size_t, fwdpp::uint_t>>;

#define GENOTYPE_MATRIX(POPTYPE, CLASSTYPE)                                   \
    m.def("genotype_matrix",                                                  \
          [](const POPTYPE &pop, const std::vector<std::size_t> &individuals, \
             const keytype &neutral_keys, const keytype &selected_keys) {     \
              return fwdpp::genotype_matrix<POPTYPE>(                         \
                  pop, individuals, neutral_keys, selected_keys);             \
          },                                                                  \
          "Generate a :class:`fwdpy11.sampling.DataMatrix` from a "           \
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
              return fwdpp::haplotype_matrix<POPTYPE>(                        \
                  pop, individuals, neutral_keys, selected_keys);             \
          },                                                                  \
          "Generate a :class:`fwdpy11.sampling.DataMatrix` from a "           \
          ":class:`" CLASSTYPE "` object.\n"                                  \
          "The DataMatrix will be encoded as haplotypes.\n\n"                 \
          ":param pop: A population object.\n"                                \
          ":param neutral_keys: The return value from "                       \
          ":func:`fwdpy11.sampling.mutation_keys`.\n\n"                       \
          ":param selected_keys: The return value from "                      \
          ":func:`fwdpy11.sampling.mutation_keys`.\n\n"                       \
          ":rtype: :class:`fwdpy11.sampling.DataMatrix` encoded as a "        \
          "haplotype matrix\n");

    MUTATION_KEYS(fwdpy11::SlocusPop, "fwdpy11.SlocusPop");
    MUTATION_KEYS(fwdpy11::MlocusPop, "fwdpy11.MlocusPop");

    GENOTYPE_MATRIX(fwdpy11::SlocusPop, "fwdpy11.SlocusPop");
    GENOTYPE_MATRIX(fwdpy11::MlocusPop, "fwdpy11.MlocusPop");

    HAPLOTYPE_MATRIX(fwdpy11::SlocusPop, "fwdpy11.SlocusPop");
    HAPLOTYPE_MATRIX(fwdpy11::MlocusPop, "fwdpy11.MlocusPop");

    m.def("matrix_to_sample",
          [](const fwdpp::data_matrix &m)

          {
              auto neutral = fwdpy11::matrix_to_sample(m.neutral);
              auto selected = fwdpy11::matrix_to_sample(m.selected);
              return py::make_tuple(std::move(neutral), std::move(selected));
          },
          R"delim(
          Convert a :class:`fwdpy11.sampling.DataMatrix` into a tuple representing the
          neutral and selected data, resepectively, as a list of (positon, string) tuples.

          See :ref:`sampling` and :ref:`datamatrix` for more details.

          .. versionadded:: 0.1.1

          .. versionchanged:: 0.2.0
                
              The return value is now a tuple (neutral, selected)

          :param m: A data matrix
          :type m: :class:`fwdpy11.sampling.DataMatrix`

          :rtype: tuple

          :return: A tuple of the (neutral, selected) mutations. Each element of the tuple holds data in list format.
              Positions are in same order 
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
          py::arg("m"));

    m.def("separate_samples_by_loci", &fwdpy11::separate_samples_by_loci,
          R"delim(
            Convert the output from :func:`fwdpy11.sampling.matrix_to_sample` into 
            separate records per locus.

			.. versionadded:: 0.1.1

			.. versionchanged:: 0.1.2
				Take a list of positions as arguments and not a population object.
		
            .. versionchanged:: 0.2.0
                Return type changed from dict to list

            :param boundaries: A list of [start,stop) tuples representing positions.
            :param sample: The return value of :func:`fwdpy11.sampling.matrix_to_sample`

            :rtype: list

            :return: The data returned follow the same structure 
				as :func:`fwdpy11.sampling.matrix_to_sample`,
                but there is one entry per locus. The data for
                consecutive loci are consecutive elements
                in the return value.
            )delim");
}
