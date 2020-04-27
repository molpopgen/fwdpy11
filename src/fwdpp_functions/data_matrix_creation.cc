#include <fwdpp/data_matrix.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/sampling/data_matrix_functions.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


void
init_data_matrix_creation_functions(py::module &m)
{
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

    MUTATION_KEYS(fwdpy11::DiploidPopulation, "fwdpy11.DiploidPopulation");

    GENOTYPE_MATRIX(fwdpy11::DiploidPopulation, "fwdpy11.DiploidPopulation");

    HAPLOTYPE_MATRIX(fwdpy11::DiploidPopulation, "fwdpy11.DiploidPopulation");

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
          )delim",
          py::arg("m"));
}

