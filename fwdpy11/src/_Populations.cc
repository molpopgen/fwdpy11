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
#include <fwdpp/sugar/sampling.hpp>
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
        .def(py::init<fwdpp::uint_t>(), "Construct with an unsigned integer "
                                        "representing the initial "
                                        "population size.")
        .def(py::init<const fwdpy11::SlocusPop::dipvector_t&,
                      const fwdpy11::SlocusPop::gcont_t&,
                      const fwdpy11::SlocusPop::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
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
            })
        .def(py::pickle(
            [](const fwdpy11::SlocusPop& pop) -> py::object {
                auto pb = py::bytes(
                    fwdpy11::serialization::serialize_details(&pop));
                return py::make_tuple(std::move(pb), pop.popdata,
                                      pop.popdata_user);
            },
            [](py::object pickled) -> fwdpy11::SlocusPop {
                try
                    {
                        auto s = pickled.cast<py::bytes>();
                        return fwdpy11::serialization::deserialize_details<
                            fwdpy11::SlocusPop>()(s, 1);
                    }
                catch (std::runtime_error& eas)
                    {
                        PyErr_Clear();
                    }
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 3)
                    {
                        throw std::runtime_error(
                            "expected tuple with 3 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto rv = fwdpy11::serialization::deserialize_details<
                    fwdpy11::SlocusPop>()(s, 1);
                rv.popdata = t[1];
                rv.popdata_user = t[2];
                return rv;
            }))
        .def("sample",
             [](const fwdpy11::SlocusPop& pop, const bool separate,
                const bool remove_fixed, py::kwargs kwargs) -> py::object {
                 py::object rv;

                 std::vector<std::size_t> ind = get_individuals(pop.N, kwargs);

                 if (separate)
                     {
                         auto temp
                             = fwdpp::sample_separate(pop, ind, remove_fixed);
                         rv = py::make_tuple(temp.first, temp.second);
                     }
                 else
                     {
                         auto temp = fwdpp::sample(pop, ind, remove_fixed);
                         py::list tlist = py::cast(temp);
                         rv = tlist;
                     }
                 return rv;
             },
             py::arg("separate") = true, py::arg("remove_fixed") = true,
             R"delim(
             Sample diploids from the population.

             :param separate: (True) Return neutral and selected variants separately.
             :param remove_fixed: (True) Remove variants fixed in the sample.
             :param kwargs: See below.

             :rtype: object

             :returns: Haplotype information for a sample.

             The valid kwargs are:

             * individuals, which should be a list of non-negative integers
             * rng, which should be a :class:`fwdpy11.GSLrng`
             * nsam, which should be a positive integer
             
             The latter two kwargs must be used together, and will generate a sample of
             ``nsam`` individuals taken *with replacement* from the population. 

             The return value is structured around a list of tuples.  Each tuple
             is (position, genotype), where genotype are encoded as 0/1 = ancestral/derived.
             From index 0 to 2*nsam - 1 (or 2*len(individuals) -1), adjacent pairs of 
             values represent diploid genotype data.  Across sites, the data represent
             individual haplotypes.

             When `separate` is `True`, a tuple of two such lists is returned.  The first
             list is for genotypes at neutral variants.  The second list is for non-neutral
             variants.

			 .. versionadded:: 0.1.4
             )delim")
        .def("add_mutations", &fwdpy11::SlocusPop::add_mutations);

    py::class_<fwdpy11::MlocusPop, fwdpy11::Population>(
        m, "_MlocusPop", "Representation of a multi-locus population")
        .def(py::init<fwdpp::uint_t, std::vector<std::pair<double, double>>>(),
             py::arg("N"), py::arg("locus_boundaries"))
        .def(py::init<const fwdpy11::MlocusPop::dipvector_t&,
                      const fwdpy11::MlocusPop::gcont_t&,
                      const fwdpy11::MlocusPop::mcont_t&,
                      std::vector<std::pair<double, double>>>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
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
            })
        .def(py::pickle(
            [](const fwdpy11::MlocusPop& pop) -> py::object {
                auto pb = py::bytes(
                    fwdpy11::serialization::serialize_details(&pop));
                return py::make_tuple(std::move(pb), pop.popdata,
                                      pop.popdata_user);
            },
            [](py::object pickled) -> fwdpy11::MlocusPop {
                try
                    {
                        auto s = pickled.cast<py::bytes>();
                        return fwdpy11::serialization::deserialize_details<
                            fwdpy11::MlocusPop>()(
                            s, 1,
                            std::vector<std::pair<double, double>>{
                                { 0., 1 } });
                    }
                catch (std::runtime_error& eas)
                    {
                        PyErr_Clear();
                    }
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 3)
                    {
                        throw std::runtime_error(
                            "expected tuple with 3 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto rv = fwdpy11::serialization::deserialize_details<
                    fwdpy11::MlocusPop>()(
                    s, 1, std::vector<std::pair<double, double>>{ { 0., 1 } });
                rv.popdata = t[1];
                rv.popdata_user = t[2];
                return rv;
            }))
        .def("sample",
             [](const fwdpy11::MlocusPop& pop, const bool separate,
                const bool remove_fixed, py::kwargs kwargs) -> py::list {
                 py::list rv;

                 std::vector<std::size_t> ind = get_individuals(pop.N, kwargs);

                 if (separate)
                     {
                         auto temp
                             = fwdpp::sample_separate(pop, ind, remove_fixed);
                         rv = py::cast(temp);
                     }
                 else
                     {
                         auto temp = fwdpp::sample(pop, ind, remove_fixed);
                         rv = py::cast(temp);
                     }
                 return rv;
             },
             R"delim(
            Sample diploids from the population.

             :param separate: (True) Return neutral and selected variants separately.
             :param remove_fixed: (True) Remove variants fixed in the sample.
             :param kwargs: See below.

             :rtype: list

             :returns: Haplotype information for a sample.

             The valid kwargs are:

             * individuals, which should be a list of non-negative integers
             * rng, which should be a :class:`fwdpy11.GSLrng`
             * nsam, which should be a positive integer
             
             The latter two kwargs must be used together, and will generate a sample of
             ``nsam`` individuals taken *with replacement* from the population. 

             The return value is structured around a list of tuples.  Each tuple
             is (position, genotype), where genotype are encoded as 0/1 = ancestral/derived.
             From index 0 to 2*nsam - 1 (or 2*len(individuals) -1), adjacent pairs of 
             values represent diploid genotype data.  Across sites, the data represent
             individual haplotypes.

             When `separate` is `True`, a list of tuples is returned.  The length of this
             list equals the number of loci, and each tuple has two elements. The first
             element is a list genotypes at neutral variants, following the format described
             above.   The second element is for non-neutral variants.

			 .. versionadded:: 0.1.4
             )delim",
             py::arg("separate") = true, py::arg("remove_fixed") = true)
        .def("add_mutations", &fwdpy11::MlocusPop::add_mutations);
}
