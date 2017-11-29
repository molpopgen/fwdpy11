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
#include <pybind11/pytypes.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types.hpp>
#include <fwdpp/fwdpp/sugar/sampling.hpp>

namespace py = pybind11;

using fwdpp_popgenmut_base = fwdpy11::singlepop_t::popbase_t;
using singlepop_sugar_base = fwdpy11::singlepop_t::base;
using multilocus_sugar_base = fwdpy11::multilocus_t::base;
using multilocus_popgenmut_base = multilocus_sugar_base::popbase_t;
using singlepop_generalmut_vec_sugar_base = fwdpy11::singlepop_gm_vec_t::base;
using singlepop_generalmut_vec_base
    = singlepop_generalmut_vec_sugar_base::popbase_t;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::dipvector_t>);
PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::uint_t>);
PYBIND11_MAKE_OPAQUE(
    std::vector<double>); // for generalmut_vec::s and generalmut_vec::h

struct flattened_popgenmut
{
    KTfwd::uint_t g;
    decltype(KTfwd::popgenmut::xtra) label;
    std::int8_t neutral;
    double pos, s, h;
};

inline flattened_popgenmut
make_flattened_popgenmut(const KTfwd::popgenmut& m)
{
    flattened_popgenmut rv;
    rv.g = m.g;
    rv.label = m.xtra;
    rv.neutral = m.neutral;
    rv.pos = m.pos;
    rv.s = m.s;
    rv.h = m.h;
    return rv;
}

struct diploid_traits
{
    double g, e, w;
};

inline diploid_traits
make_diploid_traits(const fwdpy11::diploid_t& dip)
{
    diploid_traits d;
    d.g = dip.g;
    d.e = dip.e;
    d.w = dip.w;
    return d;
}

struct diploid_gametes
{
    std::size_t locus, first, second;
};

inline diploid_gametes
make_diploid_gametes(const fwdpy11::diploid_t& dip, const std::size_t locus)
{
    diploid_gametes d;
    d.locus = locus;
    d.first = dip.first;
    d.second = dip.second;
    return d;
}

PYBIND11_MAKE_OPAQUE(std::vector<flattened_popgenmut>);
PYBIND11_MAKE_OPAQUE(std::vector<diploid_traits>);
PYBIND11_MAKE_OPAQUE(std::vector<diploid_gametes>);

namespace
{
    static const auto MCOUNTS_DOCSTRING = R"delim(
    List of number of occurrences of elements in 
    a population objecst "mutations" container.

    The values are unsigned 32-bit integers.  

    .. note::
        Some values may be 0.  These represent *extinct* variants.  You will typically want to avoid processing such mutations.
)delim";

    static const auto DIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.fwdpy11_types.DiploidContainer`.
   )delim";

    static const auto FIXATIONS_DOCSTRING
        = R"delim(A :class:`fwdpy11.fwdpp_types.MutationContainer` of fixed variants.)delim";

    static const auto FIXATION_TIMES_DOCSTRING =
        R"delim(A list of fixation times corresponding to the elements in "fixations" for this type.)delim";

    static const auto GAMETES_DOCSTRING
        = R"delim(A :class:`fwdpy11.fwdpp_types.GameteContainer`.)delim";

    static const auto MUTATIONS_DOCSTRING = R"delim(
    List of :class:`fwdpy11.fwdpp_types.Mutation`.

    .. note:: 
        This list contains **both** extinct *and* extant mutations.  
        To distinguish them, use the locations of nonzero values in "mcounts" 
        for an instance of this type."
    )delim";
}

PYBIND11_MODULE(fwdpy11_types, m)
{
    m.doc() = "Wrap C++ types specific to fwdpy11.";

    py::class_<fwdpy11::GSLrng_t>(
        m, "GSLrng", "Random number generator based on a mersenne twister.")
        .def(py::init<unsigned>(),
             "Constructor takes unsigned integer as a seed");

    py::class_<fwdpy11::diploid_t>(
        m, "SingleLocusDiploid",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_static("create", &fwdpy11::diploid_t::create)
        .def_readonly("first", &fwdpy11::diploid_t::first,
                      "Key to first gamete. (read-only)")
        .def_readonly("second", &fwdpy11::diploid_t::second,
                      "Key to second gamete. (read-only)")
        .def_readonly("w", &fwdpy11::diploid_t::w, "Fitness. (read-only)")
        .def_readonly("g", &fwdpy11::diploid_t::g,
                      "Genetic value (read-only).")
        .def_readonly("e", &fwdpy11::diploid_t::e,
                      "Random/environmental effects (read-only).")
        .def_readonly("label", &fwdpy11::diploid_t::label,
                      "Index of the diploid in its deme")
        .def_readonly("parental_data", &fwdpy11::diploid_t::parental_data,
                      R"delim(
				A tuple of the parental labels.

				.. versionadded:: 0.1.4

				.. note::
					This field is not pickled. The representation
					will change in future releases.
				)delim")
        .def(py::pickle(
            [](const fwdpy11::diploid_t& d) {
                return py::make_tuple(d.first, d.second, d.w, d.g, d.e,
                                      d.label);
            },
            [](py::tuple t) {
                std::unique_ptr<fwdpy11::diploid_t> d(new fwdpy11::diploid_t(
                    t[0].cast<std::size_t>(), t[1].cast<std::size_t>()));
                d->w = t[2].cast<double>();
                d->g = t[3].cast<double>();
                d->e = t[4].cast<double>();
                d->label = t[5].cast<decltype(fwdpy11::diploid_t::label)>();
                return d;
            }));

    py::bind_vector<fwdpy11::dipvector_t>(
        m, "DiploidContainer", py::module_local(false),
        "C++ representation of a list of "
        ":class:`fwdpy11.fwdpy11_types."
        "SingleLocusDiploid`.  Typically, access will be read-only.")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids) {
                 std::vector<diploid_traits> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.push_back(make_diploid_traits(dip));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_traits> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         // range-check here
                         auto&& dip = diploids.at(r(i));
                         rv.push_back(make_diploid_traits(dip));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids, py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_traits> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start);
                         rv.push_back(make_diploid_traits(dip));
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids) {
                 std::vector<diploid_gametes> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.push_back(make_diploid_gametes(dip, 0));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_gametes> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         rv.push_back(make_diploid_gametes(diploids.at(i), 0));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids, py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_gametes> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         rv.push_back(
                             make_diploid_gametes(diploids.at(start), 0));
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim");

    py::bind_vector<std::vector<fwdpy11::dipvector_t>>(
        m, "VecDiploidContainer", py::module_local(false),
        "Vector of "
        ":class:`fwdpy11.fwdpy11_types."
        "SingleLocusDiploid`.")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                 std::vector<diploid_traits> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.emplace_back(make_diploid_traits(dip.at(0)));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                 std::vector<diploid_gametes> rv;
                 std::size_t locus;
                 for (auto&& dip : diploids)
                     {
                         locus = 0;
                         for (auto&& di : dip)
                             {
                                 rv.emplace_back(
                                     make_diploid_gametes(di, locus++));
                             }
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_traits> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         // range-check here
                         auto&& dip = diploids.at(r(i)).at(0);
                         rv.push_back(make_diploid_traits(dip));
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_gametes> rv;
                 rv.reserve(r.shape(0));
                 std::size_t locus;
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         auto&& dip = diploids.at(r(i));
                         locus = 0;
                         for (auto&& di : dip)
                             {
                                 rv.emplace_back(
                                     make_diploid_gametes(di, locus++));
                             }
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_traits> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start).at(0);
                         rv.push_back(make_diploid_traits(dip));
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_gametes> rv;
                 rv.reserve(slicelength);
                 std::size_t locus;
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start).at(0);
                         locus = 0;
                         rv.push_back(make_diploid_gametes(dip, locus++));
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.fwdpy11_types.VecDipGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim");

    py::bind_vector<std::vector<KTfwd::uint_t>>(
        m, "VectorUint32", "Vector of unsigned 32-bit integers.",
        py::buffer_protocol(), py::module_local(false));
    py::bind_vector<fwdpy11::gcont_t>(m, "GameteContainer",
                                      py::module_local(false),
                                      "C++ representations of a list of "
                                      ":class:`fwdpy11.fwdpp_types.Gamete`.  "
                                      "Typically, access will be read-only.");

    PYBIND11_NUMPY_DTYPE(flattened_popgenmut, pos, s, h, g, label, neutral);
    PYBIND11_NUMPY_DTYPE(diploid_traits, g, e, w);
    PYBIND11_NUMPY_DTYPE(diploid_gametes, locus, first, second);
    py::bind_vector<std::vector<flattened_popgenmut>>(
        m, "VecMutStruct", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the data fields in a "
        ":class:`fwdpy11.fwdpp_types.Mutation`.

        .. versionadded: 0.1.2
        )delim");

    py::bind_vector<std::vector<diploid_traits>>(
        m, "VecDipTraits", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the g,e,w data fields in a "
        ":class:`fwdpy11.fwdpp_types.SingleLocusDiploid`.

        .. versionadded: 0.1.2
        )delim");

    py::bind_vector<std::vector<diploid_gametes>>(
        m, "VecDipGametes", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the locus and gamete index data fields in a "
        ":class:`fwdpy11.fwdpp_types.SingleLocusDiploid`.

        .. versionadded: 0.1.2
        )delim");

    py::bind_vector<fwdpy11::mcont_t>(
        m, "MutationContainer", "C++ representation of a list of "
                                ":class:`fwdpy11.fwdpp_types.Mutation`.  "
                                "Typically, access will be read-only.",
        py::module_local(false))
        .def("array",
             [](const fwdpy11::mcont_t& mc) {
                 std::vector<flattened_popgenmut> rv;
                 rv.reserve(mc.size());
                 for (auto&& m : mc)
                     {
                         flattened_popgenmut t;
                         t.pos = m.pos;
                         t.s = m.s;
                         t.h = m.h;
                         t.g = m.g;
                         t.label = m.xtra;
                         t.neutral = m.neutral;
                         rv.push_back(std::move(t));
                     }
                 return rv;
             },
             R"delim(
        :rtype: :class:`fwdpy11.fwdpy11_types.VecMutStruct`.
        
        The return value should be coerced into a Numpy 
        array for processing.

        .. versionadded: 0.1.2
        )delim");

    // expose the base classes for population types
    py::class_<fwdpp_popgenmut_base>(m, "SlocusPopMutationBase");
    py::class_<multilocus_popgenmut_base>(m, "MlocusMutationBase");
    py::class_<singlepop_sugar_base, fwdpp_popgenmut_base>(m, "SinglepopBase");
    py::class_<multilocus_sugar_base, multilocus_popgenmut_base>(m,
                                                                 "MlocusBase");

    py::class_<singlepop_generalmut_vec_base>(m, "SlocusPopGeneralMutVecBase");
    py::class_<singlepop_generalmut_vec_sugar_base,
               singlepop_generalmut_vec_base>(
        m, "SlocusPopGeneralMutVecSugarBase");

    // Expose the type based on fwdpp's "sugar"
    // layer
    py::class_<fwdpy11::singlepop_t, singlepop_sugar_base>(
        m, "SlocusPop", "Population object representing a single "
                        "deme and a "
                        "single genomic region.")
        .def(py::init<unsigned>(), "Construct with an unsigned integer "
                                   "representing the initial "
                                   "population size.")
        .def(py::init<const fwdpy11::singlepop_t::dipvector_t&,
                      const fwdpy11::singlepop_t::gcont_t&,
                      const fwdpy11::singlepop_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             ..versionadded:: 0.1.4
             )delim")
        .def_static("create", &fwdpy11::singlepop_t::create,
                    py::arg("diploids"), py::arg("gametes"),
                    py::arg("mutations"),
                    R"delim(
                    Create a new object from input data.
                    Unlike the constructor method, this method results
                    in no temporary copies of input data.
                    
                    :param diplods: A :class:`fwdpy11.fwdpy11_types.DiploidContainer`
                    :param gametes: A :class:`fwdpy11.fwdpy11_types.GameteContainer`
                    :param mutations: A :class:`fwdpy11.fwdpy11_types.MutationContainer`

                    :rtype: :class:`fwdpy11.fwdpy11_types.SlocusPop`

                    .. versionadded:: 0.1.4

                    .. note::
                        See :ref:`popobjects` for example use.
                    )delim")
        .def_static("create_with_fixations",
                    &fwdpy11::singlepop_t::create_with_fixations)
        .def("clear", &fwdpy11::singlepop_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::singlepop_t::generation,
                      R"delim(
                      The current generation. A population starts at 
                      generation 0.
                      )delim")
        .def_readonly("N", &fwdpy11::singlepop_t::N,
                      R"delim(
                      The current population size.
                      )delim")
        .def_readonly("diploids", &fwdpy11::singlepop_t::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_readonly("mutations", &fwdpp_popgenmut_base::mutations,
                      MUTATIONS_DOCSTRING)
        .def_readonly("mcounts", &fwdpp_popgenmut_base::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly("fixations", &fwdpp_popgenmut_base::fixations,
                      FIXATIONS_DOCSTRING)
        .def_readonly("fixation_times", &fwdpp_popgenmut_base::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def_readonly("gametes", &fwdpp_popgenmut_base::gametes,
                      GAMETES_DOCSTRING)
        .def(py::pickle(
            [](const fwdpy11::singlepop_t& pop) {
                return py::bytes(pop.serialize());
            },
            [](py::bytes s) {
                return std::unique_ptr<fwdpy11::singlepop_t>(
                    new fwdpy11::singlepop_t(s));
            }))
        .def("__eq__",
             [](const fwdpy11::singlepop_t& lhs,
                const fwdpy11::singlepop_t& rhs) { return lhs == rhs; })
        .def("sample",
             [](const fwdpy11::singlepop_t& pop, const fwdpy11::GSLrng_t& rng,
                const std::int64_t nsam, const bool separate,
                const bool remove_fixed) -> py::object {
                 if (nsam <= 0)
                     {
                         throw std::invalid_argument(
                             "sample size must be > 0");
                     }
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(rng.get(), pop, nsam,
                                                         remove_fixed);
                         auto t = py::make_tuple(std::move(s.first),
                                                 std::move(s.second));
                         return t;
                     }
                 py::list rv;
                 auto s = KTfwd::sample(rng.get(), pop, nsam, remove_fixed);
                 for (auto& i : s)
                     {
                         rv.append(i);
                     }
                 return rv;
             },
             py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param nsam: An integer representing the sample size in chromosomes (2 times number of diploids).
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.

             .. note::
                If you want sampling *without* replacement, see
                :func:`~fwdpy11.fwdpy11_types.SlocusPop.sample_ind`.
             )delim")
        .def("sample_ind",
             [](const fwdpy11::singlepop_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(pop, individuals,
                                                         remove_fixed);
                         auto t = py::make_tuple(std::move(s.first),
                                                 std::move(s.second));
                         return t;
                     }
                 py::list rv;
                 auto s = KTfwd::sample(pop, individuals, remove_fixed);
                 for (auto& i : s)
                     {
                         rv.append(i);
                     }
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param individuals: The individuals to sample.
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The final sample size is ``2*len(individuals)``.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.
             )delim");

    py::class_<fwdpy11::multilocus_t, multilocus_sugar_base>(
        m, "MlocusPop", "Representation of a multi-locus, single "
                        "deme system.")
        .def(py::init<unsigned, unsigned>(), py::arg("N"), py::arg("nloci"),
             "Construct with population size and "
             "number of loci.")
        .def(py::init<const fwdpy11::multilocus_t::dipvector_t&,
                      const fwdpy11::multilocus_t::gcont_t&,
                      const fwdpy11::multilocus_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             ..versionadded:: 0.1.4
             )delim")
        .def(py::init<unsigned, unsigned,
                      const std::vector<std::pair<double, double>>&>(),
             py::arg("N"), py::arg("nloci"), py::arg("locus_boundaries"))
        .def_static("create", &fwdpy11::multilocus_t::create,
                    py::arg("diploids"), py::arg("gametes"),
                    py::arg("mutations"),
                    R"delim(
                    Create a new object from input data.
                    Unlike the constructor method, this method results
                    in no temporary copies of input data.
                    
                    :param diplods: A :class:`fwdpy11.fwdpy11_types.VecDiploidContainer`
                    :param gametes: A :class:`fwdpy11.fwdpy11_types.GameteContainer`
                    :param mutations: A :class:`fwdpy11.fwdpy11_types.MutationContainer`

                    :rtype: :class:`fwdpy11.fwdpy11_types.MlocusPop`

                    .. versionadded:: 0.1.4

                    .. note::
                        See :ref:`popobjects` for example use.
                    )delim")
        .def_static("create_with_fixations",
                    &fwdpy11::multilocus_t::create_with_fixations)
        .def("clear", &fwdpy11::multilocus_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::multilocus_t::generation,
                      "The current generation.")
        .def_readonly("N", &fwdpy11::multilocus_t::N,
                      "Curent population size.")
        .def_readonly("nloci", &fwdpy11::multilocus_t::nloci,
                      "Number of loci.")
        .def_readonly("diploids", &fwdpy11::multilocus_t::diploids)
        .def_readonly("mutations", &fwdpy11::multilocus_t::mutations,
                      MUTATIONS_DOCSTRING)
        .def_readonly("gametes", &fwdpy11::multilocus_t::gametes,
                      GAMETES_DOCSTRING)
        .def_readonly("mcounts", &fwdpy11::multilocus_t::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly("fixations", &fwdpy11::multilocus_t::fixations,
                      FIXATIONS_DOCSTRING)
        .def_readonly("fixation_times", &fwdpy11::multilocus_t::fixation_times,
                      FIXATIONS_DOCSTRING)
        .def_readwrite("locus_boundaries",
                       &fwdpy11::multilocus_t::locus_boundaries,
                       "[beg,end) positions for each locus")
        .def(py::pickle(
            [](const fwdpy11::multilocus_t& pop) -> py::object {
                auto pb = py::bytes(pop.serialize());
                py::list pdata;
                for (auto& d : pop.diploids)
                    {
                        pdata.append(d[0].parental_data);
                    }
                return py::make_tuple(std::move(pb), std::move(pdata));
            },
            [](py::object pickled) {
                
			
			try
                    {
                        auto s = pickled.cast<py::bytes>();
                        return std::unique_ptr<fwdpy11::multilocus_t>(
                            new fwdpy11::multilocus_t(s));
                    }
                catch (std::runtime_error& eas)
                    {
                        try
                            {
                                auto t = pickled.cast<py::tuple>();
                                auto s = t[0].cast<py::bytes>();
                                auto l = t[1].cast<py::list>();
                                auto rv
                                    = std::unique_ptr<fwdpy11::multilocus_t>(
                                        new fwdpy11::multilocus_t(s));
                                for (std::size_t i = 0;
                                     i < rv->diploids.size(); ++i)
                                    {
                                        rv->diploids[i][0].parental_data
                                            = l[i];
                                    }
                                return rv;
                            }
                        catch (py::error_already_set& eas)
                            {
								eas.restore();
                                auto t = pickled.cast<py::tuple>();
                                auto s = t[0].cast<py::bytes>();
                                auto l = t[1].cast<py::list>();
                                auto rv
                                    = std::unique_ptr<fwdpy11::multilocus_t>(
                                        new fwdpy11::multilocus_t(s));
                                for (std::size_t i = 0;
                                     i < rv->diploids.size(); ++i)
                                    {
                                        rv->diploids[i][0].parental_data
                                            = l[i];
                                    }
                                return rv;
                            }
                    }
                                throw std::runtime_error("bad juju!");
            }))
        .def("__eq__",
             [](const fwdpy11::multilocus_t& lhs,
                const fwdpy11::multilocus_t& rhs) { return lhs == rhs; })
        .def(
            "sample",
            [](const fwdpy11::multilocus_t& pop, const fwdpy11::GSLrng_t& rng,
               const std::int64_t nsam, const bool separate,
               const bool remove_fixed) -> py::object {
                if (nsam <= 0)
                    {
                        throw std::invalid_argument("sample size must be > 0");
                    }
                if (separate)
                    {
                        auto s = KTfwd::sample_separate(rng.get(), pop, nsam,
                                                        remove_fixed);
                        py::list rv;
                        for (auto& i : s)
                            {
                                rv.append(py::make_tuple(std::move(i.first),
                                                         std::move(i.second)));
                            }
                        return rv;
                    }
                py::list rv;
                auto s = KTfwd::sample(rng.get(), pop, nsam, remove_fixed);
                for (auto& i : s)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
            py::arg("remove_fixed") = true,
            R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param nsam: An integer representing the sample size in chromosomes (2 times number of diploids).
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.

             .. note::
                If you want sampling *without* replacement, see
                :func:`~fwdpy11.fwdpy11_types.MlocusPop.sample_ind`.
             )delim")
        .def("sample_ind",
             [](const fwdpy11::multilocus_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(pop, individuals,
                                                         remove_fixed);
                         py::list rv;
                         for (auto& i : s)
                             {
                                 rv.append(py::make_tuple(
                                     std::move(i.first), std::move(i.second)));
                                 return rv;
                             }
                     }
                 py::list rv;
                 auto s = KTfwd::sample(pop, individuals, remove_fixed);
                 for (auto& i : s)
                     {
                         rv.append(i);
                     }
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param individuals: The individuals to sample.
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The final sample size is ``2*len(individuals)``.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.
             )delim");

    py::class_<fwdpy11::singlepop_gm_vec_t,
               singlepop_generalmut_vec_sugar_base>(
        m, "SlocusPopGeneralMutVec",
        "Single-deme object using "
        ":class:`fwpy11.fwdpp_types.GeneralMutVec`"
        " as "
        "the mutation type.")
        .def(py::init<unsigned>(), py::arg("N"),
             "Construct object with N diploids.")
        .def(py::init<const fwdpy11::singlepop_gm_vec_t::dipvector_t&,
                      const fwdpy11::singlepop_gm_vec_t::gcont_t&,
                      const fwdpy11::singlepop_gm_vec_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             ..versionadded:: 0.1.4
             )delim")
        .def_static("create", &fwdpy11::singlepop_gm_vec_t::create,
                    py::arg("diploids"), py::arg("gametes"),
                    py::arg("mutations"),
                    R"delim(
                    Create a new object from input data.
                    Unlike the constructor method, this method results
                    in no temporary copies of input data.
                    
                    :param diplods: A :class:`fwdpy11.fwdpy11_types.DiploidContainer`
                    :param gametes: A :class:`fwdpy11.fwdpy11_types.GameteContainer`
                    :param mutations: A :class:`fwdpy11.fwdpp_types.VectorGeneralMutVec`

                    :rtype: :class:`fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec`

                    .. versionadded:: 0.1.4

                    .. note::
                        See :ref:`popobjects` for example use.
                    )delim")
        .def_static("create_with_fixations",
                    &fwdpy11::singlepop_gm_vec_t::create_with_fixations)
        .def("clear", &fwdpy11::singlepop_gm_vec_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::singlepop_gm_vec_t::generation,
                      "The current generation.")
        .def_readonly("N", &fwdpy11::singlepop_gm_vec_t::N,
                      "Curent population size.")
        .def_readonly("diploids", &fwdpy11::singlepop_gm_vec_t::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_readonly("mutations", &fwdpy11::singlepop_gm_vec_t::mutations,
                      "A list of "
                      ":class:`fwdpy11.fwdpp_types."
                      "VectorGeneralMutVec`.")
        .def_readonly("gametes", &fwdpy11::singlepop_gm_vec_t::gametes,
                      GAMETES_DOCSTRING)
        .def_readonly("mcounts", &fwdpy11::singlepop_gm_vec_t::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly("fixations", &fwdpy11::singlepop_gm_vec_t::fixations,
                      "A list of "
                      ":class:`fwdpy11.fwdpp_types."
                      "VectorGeneralMutVec`.")
        .def_readonly("fixation_times",
                      &fwdpy11::singlepop_gm_vec_t::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def(py::pickle(
            [](const fwdpy11::singlepop_gm_vec_t& pop) {
                return py::bytes(pop.serialize());
            },
            [](py::bytes s) {
                return std::unique_ptr<fwdpy11::singlepop_gm_vec_t>(
                    new fwdpy11::singlepop_gm_vec_t(s));
            }))
        .def("__eq__",
             [](const fwdpy11::singlepop_gm_vec_t& lhs,
                const fwdpy11::singlepop_gm_vec_t& rhs) { return lhs == rhs; })
        .def("sample",
             [](const fwdpy11::singlepop_gm_vec_t& pop,
                const fwdpy11::GSLrng_t& rng, const std::int64_t nsam,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (nsam <= 0)
                     {
                         throw std::invalid_argument(
                             "sample size must be > 0");
                     }
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(rng.get(), pop, nsam,
                                                         remove_fixed);
                         auto t = py::make_tuple(std::move(s.first),
                                                 std::move(s.second));
                         return t;
                     }
                 py::list rv;
                 auto s = KTfwd::sample(rng.get(), pop, nsam, remove_fixed);
                 for (auto& i : s)
                     {
                         rv.append(i);
                     }
                 return rv;
             },
             py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param nsam: An integer representing the sample size in chromosomes (2 times number of diploids).
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.

             .. note::
                If you want sampling *without* replacement, see
                :func:`~fwdpy11.fwdpy11_types.SlocusPopGeneralMutVec.sample_ind`.
             )delim")
        .def("sample_ind",
             [](const fwdpy11::singlepop_gm_vec_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(pop, individuals,
                                                         remove_fixed);
                         auto t = py::make_tuple(std::move(s.first),
                                                 std::move(s.second));
                         return t;
                     }
                 py::list rv;
                 auto s = KTfwd::sample(pop, individuals, remove_fixed);
                 for (auto& i : s)
                     {
                         rv.append(i);
                     }
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
             :param individuals: The individuals to sample.
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The final sample size is ``2*len(individuals)``.

             The output is a list of tuples.  Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a tuple of two lists.
             The first list is for neutral variants, and the second for non-neutral.
             )delim");
}
