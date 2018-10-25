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

// TODO: make a versionchanged entry for all things affected by "length"

#include <limits>
#include <type_traits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/types/create_pops.hpp>
#include <fwdpy11/serialization.hpp>
#include <fwdpy11/serialization/Mutation.hpp>
#include <fwdpy11/serialization/Diploid.hpp>
#include "get_individuals.hpp"

static_assert(std::is_same<fwdpy11::SlocusPop::gcont_t,
                           fwdpy11::MlocusPop::gcont_t>::value,
              "SlocusPop and MlocusPop must have same gamete container type");

static_assert(
    std::is_same<fwdpy11::SlocusPop::mcont_t,
                 fwdpy11::MlocusPop::mcont_t>::value,
    "SlocusPop and MlocusPop must have same mutation container type");

namespace py = pybind11;

namespace
{
    static const auto DIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.VecDiploid`.
   )delim";

    static const auto MLDIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.VecVecDiploid`.
   )delim";
} // namespace

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidGenotype>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<fwdpy11::DiploidGenotype>>);
PYBIND11_MAKE_OPAQUE(fwdpy11::SlocusPop::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::SlocusPop::mcont_t);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::uint_t>);

PYBIND11_MODULE(_Populations, m)
{

    py::object base_class_module
        = (pybind11::object)pybind11::module::import("fwdpy11._Population");

    py::object data_matrix_python_representation
        = (pybind11::object)py::module::import("fwdpy11.sampling")
              .attr("DataMatrix");

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

    py::class_<fwdpy11::SlocusPop, fwdpy11::Population>(
        m, "_SlocusPop", "Representation of a single-locus population")
        .def(py::init<fwdpp::uint_t, double>(),
             "Construct with an unsigned integer "
             "representing the initial "
             "population size.",
             py::arg("N"),
             py::arg("length") = std::numeric_limits<double>::max())
        .def(py::init<const fwdpy11::SlocusPop::dipvector_t&,
                      const fwdpy11::SlocusPop::gcont_t&,
                      const fwdpy11::SlocusPop::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim",
             py::arg("diploids"), py::arg("gametes"), py::arg("mutations"))
        .def(py::init<const fwdpy11::SlocusPop&>(),
             R"delim(
                Copy constructor

                .. versionadded:: 0.1.4
                )delim")
        .def("clear", &fwdpy11::SlocusPop::clear,
             "Clears all population data.")
        .def("__eq__",
             [](const fwdpy11::SlocusPop& lhs, const fwdpy11::SlocusPop& rhs) {
                 return lhs == rhs;
             })
        .def_readonly("diploids", &fwdpy11::SlocusPop::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_static(
            "create",
            [](fwdpy11::SlocusPop::dipvector_t& diploids,
               fwdpy11::SlocusPop::gcont_t& gametes,
               fwdpy11::SlocusPop::mcont_t& mutations,
               py::tuple args) -> fwdpy11::SlocusPop {
                if (args.size() == 0)
                    {
                        return fwdpy11::create_wrapper<fwdpy11::SlocusPop>()(
                            std::move(diploids), std::move(gametes),
                            std::move(mutations));
                    }
                auto& fixations = args[0].cast<fwdpy11::SlocusPop::mcont_t&>();
                auto& ftimes = args[1].cast<std::vector<fwdpp::uint_t>&>();
                auto g = args[2].cast<fwdpp::uint_t>();
                return fwdpy11::create_wrapper<fwdpy11::SlocusPop>()(
                    std::move(diploids), std::move(gametes),
                    std::move(mutations), std::move(fixations),
                    std::move(ftimes), g);
            },
            py::arg("diploids"), py::arg("gametes"), py::arg("mutations"),
            py::arg("args"))
        .def(py::pickle(
            [](const fwdpy11::SlocusPop& pop) -> py::object {
                auto pb = py::bytes(
                    fwdpy11::serialization::serialize_details(&pop));
                return pb;
            },
            [](py::object pickled) -> fwdpy11::SlocusPop {
                auto s = pickled.cast<py::bytes>();
                return fwdpy11::serialization::deserialize_details<
                    fwdpy11::SlocusPop>()(s, 1,
                                          std::numeric_limits<double>::max());
            }))
        .def("sample",
             [](const fwdpy11::SlocusPop& pop,
                const std::vector<std::size_t>& individuals,
                const bool haplotype, const bool remove_fixed) {
                 return pop.sample_individuals(individuals, haplotype,
                                               remove_fixed);
             },
             R"delim(
             Return a sample from a population.

             :param individuals: Indexes of individuals in the sample
             :type individuals: list
             :param haplotype: (True) Determines the encoding of the return type
             :type haplotype: bool
             :param remove_fixed: (True) Remove fixed variants from the sample
             :type remove_fixed: bool

             :returns: A haplotype matrix if haplotype is True. Otherwise, a genotype matrix.
             :rtype: :class:`fwdpy11.sampling.DataMatrix`

             .. versionadded:: 0.1.4

             .. versionchanged:: 0.2.0

                Change to an overloaded function returning a 
                :class:`fwdpy11.sampling.DataMatrix` instead of 
                the "classic libsequence" layout.
             )delim",
             py::arg("individuals"), py::arg("haplotype") = true,
             py::arg("remove_fixed") = true)
        .def("sample",
             [](const fwdpy11::SlocusPop& pop, const fwdpy11::GSLrng_t& rng,
                const std::uint32_t nsam, const bool haplotype,
                const bool remove_fixed) {
                 return pop.sample_random_individuals(rng, nsam, haplotype,
                                                      remove_fixed);
             },
             py::arg("rng"), py::arg("nsam"), py::arg("haplotype") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Return a sample from a population.

             :param rng: Random number generator
             :type rng: :class:`fwdpy11.GSLrng`
             :param nsam: The sample size
             :type nsam: int
             :param haplotype: (True) Determines the encoding of the return type
             :type haplotype: bool
             :param remove_fixed: (True) Remove fixed variants from the sample
             :type remove_fixed: bool

             :returns: A haplotype matrix if haplotype is True. Otherwise, a genotype matrix.
             :rtype: :class:`fwdpy11.sampling.DataMatrix`

             .. versionadded:: 0.1.4

             .. versionchanged:: 0.2.0

                Change to an overloaded function returning a 
                :class:`fwdpy11.sampling.DataMatrix` instead of 
                the "classic libsequence" layout.
             )delim")

        .def("add_mutations", &fwdpy11::SlocusPop::add_mutations);

    py::class_<fwdpy11::MlocusPop, fwdpy11::Population>(
        m, "_MlocusPop", "Representation of a multi-locus population")
        .def(py::init<fwdpp::uint_t, std::vector<std::pair<double, double>>,
                      double>(),
             py::arg("N"), py::arg("locus_boundaries"),
             py::arg("length") = std::numeric_limits<double>::max())
        .def(py::init<const fwdpy11::MlocusPop::dipvector_t&,
                      const fwdpy11::MlocusPop::gcont_t&,
                      const fwdpy11::MlocusPop::mcont_t&,
                      std::vector<std::pair<double, double>>>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim",
             py::arg("diploids"), py::arg("gametes"), py::arg("mutations"),
             py::arg("locus_boundaries"))
        .def(py::init<const fwdpy11::MlocusPop&>(),
             R"delim(
                Copy constructor.

                .. versionadded:: 0.1.4
                )delim")
        .def("clear", &fwdpy11::MlocusPop::clear,
             "Clears all population data.")
        .def("__eq__",
             [](const fwdpy11::MlocusPop& lhs, const fwdpy11::MlocusPop& rhs) {
                 return lhs == rhs;
             })
        .def_readonly("nloci", &fwdpy11::MlocusPop::nloci)
        .def_property(
            "locus_boundaries",
            [](const fwdpy11::MlocusPop& pop) { return pop.locus_boundaries; },
            [](fwdpy11::MlocusPop& pop,
               std::vector<std::pair<double, double>> locus_boundaries) {
                if (locus_boundaries.size() != pop.nloci)
                    {
                        throw std::invalid_argument(
                            "incorrect number of locus boundaries");
                    }
                pop.validate_locus_boundaries(locus_boundaries);
                pop.locus_boundaries.swap(locus_boundaries);
            },
            "List of [start,stop) positions of each locus.")
        .def_readonly("diploids", &fwdpy11::MlocusPop::diploids,
                      MLDIPLOIDS_DOCSTRING)
        .def_static(
            "create",
            [](fwdpy11::MlocusPop::dipvector_t& diploids,
               fwdpy11::MlocusPop::gcont_t& gametes,
               fwdpy11::MlocusPop::mcont_t& mutations,
               std::vector<std::pair<double, double>> locus_boundaries,
               py::tuple args) -> fwdpy11::MlocusPop {
                if (args.size() == 0)
                    {
                        return fwdpy11::create_wrapper<fwdpy11::MlocusPop>()(
                            std::move(diploids), std::move(gametes),
                            std::move(mutations), std::move(locus_boundaries));
                    }
                auto& fixations = args[0].cast<fwdpy11::MlocusPop::mcont_t&>();
                auto& ftimes = args[1].cast<std::vector<fwdpp::uint_t>&>();
                auto g = args[2].cast<fwdpp::uint_t>();
                return fwdpy11::create_wrapper<fwdpy11::MlocusPop>()(
                    std::move(diploids), std::move(gametes),
                    std::move(mutations), std::move(locus_boundaries),
                    std::move(fixations), std::move(ftimes), g);
            },
            py::arg("diploids"), py::arg("gametes"), py::arg("mutations"),
            py::arg("locus_boundaries"), py::arg("args"))
        .def(py::pickle(
            [](const fwdpy11::MlocusPop& pop) -> py::object {
                auto pb = py::bytes(
                    fwdpy11::serialization::serialize_details(&pop));
                return pb;
            },
            [](py::object pickled) -> fwdpy11::MlocusPop {
                auto s = pickled.cast<py::bytes>();
                return fwdpy11::serialization::deserialize_details<
                    fwdpy11::MlocusPop>()(
                    s, 1, std::vector<std::pair<double, double>>{ { 0., 1 } },
                    std::numeric_limits<double>::max());
            }))
        .def("sample_by_locus",
             [](const fwdpy11::MlocusPop& pop,
                const std::vector<std::size_t>& individuals,
                const bool remove_fixed) {
                 return fwdpp::sample_individuals_by_window(
                     pop, individuals, pop.locus_boundaries, true, true,
                     remove_fixed);
             },
             R"delim(
             Return one data matrix per locus.
             
             :param individuals: Indexes of the individuals in the sample
             :type individuals: list
             :param remove_fixed: (True) Remove fixed variants from the sample
             :type remove_fixed: bool

             :returns: A list of haplotype matrix
             :rtype: list

             This function differs from :func:`fwdpy11.MlocusPop.sample`
             in that a list of :class:`fwdpy11.sampling.DataMatrix` are returned.
             The positional boundaries of each matrix are given by
             :attr:`fwdpy11.MlocusPop.locus_boundaries`.

             .. versionadded:: 0.2.0
             )delim",
             py::arg("individuals"), py::arg("remove_fixed") = true)
        .def("sample",
             [](const fwdpy11::MlocusPop& pop,
                const std::vector<std::size_t>& individuals,
                const bool haplotype, const bool remove_fixed) {
                 return pop.sample_individuals(individuals, haplotype,
                                               remove_fixed);
             },
             R"delim(
             Return a sample from a population.

             :param individuals: Indexes of individuals in the sample
             :type individuals: list
             :param haplotype: (True) Determines the encoding of the return type
             :type haplotype: bool
             :param remove_fixed: (True) Remove fixed variants from the sample
             :type remove_fixed: bool

             :returns: A haplotype matrix if haplotype is True. Otherwise, a genotype matrix.
             :rtype: :class:`fwdpy11.sampling.DataMatrix`

             .. versionadded:: 0.1.4

             .. versionchanged:: 0.2.0

                Change to an overloaded function returning a 
                :class:`fwdpy11.sampling.DataMatrix` instead of 
                the "classic libsequence" layout.
             )delim",
             py::arg("individuals"), py::arg("haplotype") = true,
             py::arg("remove_fixed") = true)
        .def("sample",
             [](const fwdpy11::MlocusPop& pop, const fwdpy11::GSLrng_t& rng,
                const std::uint32_t nsam, const bool haplotype,
                const bool remove_fixed) {
                 return pop.sample_random_individuals(rng, nsam, haplotype,
                                                      remove_fixed);
             },
             py::arg("rng"), py::arg("nsam"), py::arg("haplotype") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Return a sample from a population.

             :param rng: Random number generator
             :type rng: :class:`fwdpy11.GSLrng`
             :param nsam: The sample size
             :type nsam: int
             :param haplotype: (True) Determines the encoding of the return type
             :type haplotype: bool
             :param remove_fixed: (True) Remove fixed variants from the sample
             :type remove_fixed: bool

             :returns: A haplotype matrix if haplotype is True. Otherwise, a genotype matrix.
             :rtype: :class:`fwdpy11.sampling.DataMatrix`

             .. versionadded:: 0.1.4

             .. versionchanged:: 0.2.0

                Change to an overloaded function returning a 
                :class:`fwdpy11.sampling.DataMatrix` instead of 
                the "classic libsequence" layout.
             )delim")
        .def("add_mutations", &fwdpy11::MlocusPop::add_mutations);
}
