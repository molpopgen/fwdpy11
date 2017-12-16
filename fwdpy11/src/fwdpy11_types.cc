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

PYBIND11_MAKE_OPAQUE(fwdpy11::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::mcont_t);
PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::generalmut_vec>);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::diploid_t>);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::dipvector_t>);
PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::uint_t>);

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
   A :class:`fwdpy11.VecDiploid`.
   )delim";

    static const auto FIXATIONS_DOCSTRING
        = R"delim(A :class:`fwdpy11.VecMutation` of fixed variants.)delim";

    static const auto FIXATION_TIMES_DOCSTRING =
        R"delim(A list of fixation times corresponding to the elements in "fixations" for this type.)delim";

    static const auto GAMETES_DOCSTRING
        = R"delim(A :class:`fwdpy11.VecGamete`.)delim";

    static const auto MUTATIONS_DOCSTRING = R"delim(
    List of :class:`fwdpy11.Mutation`.

    .. note:: 
        This list contains **both** extinct *and* extant mutations.  
        To distinguish them, use the locations of nonzero values in "mcounts" 
        for an instance of this type."
    )delim";

    static const auto POPDATA_DOCSTRING
        = "Python object that may be written to by a simulation. Any data "
          "written should be documented by the simulation function.\n\n.. "
          "versionadded:: 0.1.4";

    static const auto POPDATA_USER_DOCSTRING
        = "A Python object with read-write access.\n\n.. versionadded:: 0.1.4";
}

PYBIND11_MODULE(fwdpy11_types, m)
{
    m.doc() = "Wrap C++ types specific to fwdpy11.";

    py::bind_vector<std::vector<KTfwd::uint_t>>(
        m, "VecUint32", "Vector of unsigned 32-bit integers.",
        py::buffer_protocol())
        .def(py::pickle(
            [](const std::vector<KTfwd::uint_t>& v) {
                py::list rv;
                for (auto&& i : v)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<KTfwd::uint_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<KTfwd::uint_t>());
                    }
                return rv;
            }));

    py::class_<fwdpy11::diploid_t>(m, "SingleLocusDiploid")
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
				Python object representing information about parents.
				The details are simulation-dependent.

				.. versionadded:: 0.1.4
				)delim")
        .def(py::pickle(
            [](const fwdpy11::diploid_t& d) {
                return py::make_tuple(d.first, d.second, d.w, d.g, d.e,
                                      d.label, d.parental_data);
            },
            [](py::tuple t) {
                std::unique_ptr<fwdpy11::diploid_t> d(new fwdpy11::diploid_t(
                    t[0].cast<std::size_t>(), t[1].cast<std::size_t>()));
                d->w = t[2].cast<double>();
                d->g = t[3].cast<double>();
                d->e = t[4].cast<double>();
                d->label = t[5].cast<decltype(fwdpy11::diploid_t::label)>();
                // Unpickle the Python object.
                // The if statement is for backwards compatibility.
                if (t.size() == 7)
                    {
                        d->parental_data = t[6];
                    }
                return d;
            }))
        .def("__eq__", [](const fwdpy11::diploid_t& a,
                          const fwdpy11::diploid_t& b) { return a == b; });

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
    py::class_<fwdpy11::singlepop_t, singlepop_sugar_base>(m, "SlocusPop")
        .def(py::init<unsigned>(), "Construct with an unsigned integer "
                                   "representing the initial "
                                   "population size.")
        .def(py::init<const fwdpy11::singlepop_t::dipvector_t&,
                      const fwdpy11::singlepop_t::gcont_t&,
                      const fwdpy11::singlepop_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
        .def(py::init<const fwdpy11::singlepop_t&>(),
             R"delim(
                Copy constructor

                .. versionadded:: 0.1.4
                )delim")
        .def_static(
            "create",
            [](fwdpy11::dipvector_t& diploids, fwdpy11::gcont_t& gametes,
               fwdpy11::mcont_t& mutations, py::tuple args) {
                if (args.size() == 0)
                    {
                        auto rv = fwdpy11::singlepop_t::create(
                            diploids, gametes, mutations);
                        return rv;
                    }
                auto& fixations = args[0].cast<fwdpy11::mcont_t&>();
                auto& ftimes = args[1].cast<std::vector<KTfwd::uint_t>&>();
                auto g = args[2].cast<KTfwd::uint_t>();
                return fwdpy11::singlepop_t::create_with_fixations(
                    diploids, gametes, mutations, fixations, ftimes, g);
            })
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
        .def_readonly("popdata", &fwdpy11::singlepop_t::popdata,
                      POPDATA_DOCSTRING)
        .def_readwrite("popdata_user", &fwdpy11::singlepop_t::popdata_user,
                       POPDATA_USER_DOCSTRING)
        .def(py::pickle(
            [](const fwdpy11::singlepop_t& pop) -> py::object {
                auto pb = py::bytes(pop.serialize());
                py::list pdata;
                for (auto& d : pop.diploids)
                    {
                        pdata.append(d.parental_data);
                    }
                return py::make_tuple(std::move(pb), std::move(pdata),
                                      pop.popdata, pop.popdata_user);
            },
            [](py::object pickled) {
                try
                    {
                        auto s = pickled.cast<py::bytes>();
                        return std::unique_ptr<fwdpy11::singlepop_t>(
                            new fwdpy11::singlepop_t(s));
                    }
                catch (std::runtime_error& eas)
                    {
                        PyErr_Clear();
                    }
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 4)
                    {
                        throw std::runtime_error(
                            "expected tuple with 4 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto l = t[1].cast<py::list>();
                auto rv = std::unique_ptr<fwdpy11::singlepop_t>(
                    new fwdpy11::singlepop_t(s));
                for (std::size_t i = 0; i < rv->diploids.size(); ++i)
                    {
                        rv->diploids[i].parental_data = l[i];
                    }
                rv->popdata = t[2];
                rv->popdata_user = t[3];
                return rv;
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
                         return py::make_tuple(s.first, s.second);
                     }
                 py::list rv = py::cast(
                     KTfwd::sample(rng.get(), pop, nsam, remove_fixed));
                 return rv;
             },
             py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.GSLrng`
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
                :func:`~fwdpy11.SlocusPop.sample_ind`.

			 .. versionadded:: 0.1.4
             )delim")
        .def("sample_ind",
             [](const fwdpy11::singlepop_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(pop, individuals,
                                                         remove_fixed);
                         return py::make_tuple(s.first, s.second);
                     }
                 py::list rv
                     = py::cast(KTfwd::sample(pop, individuals, remove_fixed));
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.GSLrng`
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

			 .. versionadded:: 0.1.4
             )delim");

    py::class_<fwdpy11::multilocus_t, multilocus_sugar_base>(m, "MlocusPop")
        .def(py::init<unsigned, unsigned>(), py::arg("N"), py::arg("nloci"),
             "Construct with population size and "
             "number of loci.")
        .def(py::init<const fwdpy11::multilocus_t::dipvector_t&,
                      const fwdpy11::multilocus_t::gcont_t&,
                      const fwdpy11::multilocus_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
        .def(py::init<const fwdpy11::multilocus_t&>(),
             R"delim(
                Copy constructor.

                .. versionadded:: 0.1.4
                )delim")
        .def(py::init<unsigned, unsigned,
                      const std::vector<std::pair<double, double>>&>(),
             py::arg("N"), py::arg("nloci"), py::arg("locus_boundaries"))
        .def_static(
            "create",
            [](std::vector<fwdpy11::multilocus_diploid_t>& diploids,
               fwdpy11::gcont_t& gametes, fwdpy11::mcont_t& mutations,
               py::tuple args) {
                if (args.size() == 0)
                    {
                        return fwdpy11::multilocus_t::create(diploids, gametes,
                                                             mutations);
                    }
                auto& fixations = args[0].cast<fwdpy11::mcont_t&>();
                auto& ftimes = args[1].cast<std::vector<KTfwd::uint_t>&>();
                auto g = args[2].cast<KTfwd::uint_t>();
                return fwdpy11::multilocus_t::create_with_fixations(
                    diploids, gametes, mutations, fixations, ftimes, g);
            })
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
        .def_readonly("popdata", &fwdpy11::multilocus_t::popdata,
                      POPDATA_DOCSTRING)
        .def_readwrite("popdata_user", &fwdpy11::multilocus_t::popdata_user,
                       POPDATA_USER_DOCSTRING)
        .def(py::pickle(
            [](const fwdpy11::multilocus_t& pop) -> py::object {
                auto pb = py::bytes(pop.serialize());
                py::list pdata;
                for (auto& d : pop.diploids)
                    {
                        pdata.append(d[0].parental_data);
                    }
                return py::make_tuple(std::move(pb), std::move(pdata),
                                      pop.popdata, pop.popdata_user);
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
                        PyErr_Clear();
                    }
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 4)
                    {
                        throw std::runtime_error(
                            "expected tuple with 4 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto l = t[1].cast<py::list>();
                auto rv = std::unique_ptr<fwdpy11::multilocus_t>(
                    new fwdpy11::multilocus_t(s));
                for (std::size_t i = 0; i < rv->diploids.size(); ++i)
                    {
                        rv->diploids[i][0].parental_data = l[i];
                    }
                rv->popdata = t[2];
                rv->popdata_user = t[3];
                return rv;
            }))
        .def("__eq__",
             [](const fwdpy11::multilocus_t& lhs,
                const fwdpy11::multilocus_t& rhs) { return lhs == rhs; })
        .def("sample",
             [](const fwdpy11::multilocus_t& pop, const fwdpy11::GSLrng_t& rng,
                const std::int64_t nsam, const bool separate,
                const bool remove_fixed) -> py::list {
                 if (nsam <= 0)
                     {
                         throw std::invalid_argument(
                             "sample size must be > 0");
                     }
                 py::list rv;
                 if (separate)
                     {
                         rv = py::cast(KTfwd::sample_separate(
                             rng.get(), pop, nsam, remove_fixed));
                         return rv;
                     }
                 rv = py::cast(
                     KTfwd::sample(rng.get(), pop, nsam, remove_fixed));
                 return rv;
             },
             py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.GSLrng`
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
                :func:`~fwdpy11.MlocusPop.sample_ind`.

			 .. versionadded:: 0.1.4
             )delim")
        .def("sample_ind",
             [](const fwdpy11::multilocus_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::list {
                 py::list rv;
                 if (separate)
                     {
                         rv = py::cast(KTfwd::sample_separate(pop, individuals,
                                                              remove_fixed));
                     }
                 rv = py::cast(KTfwd::sample(pop, individuals, remove_fixed));
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.GSLrng`
             :param individuals: The individuals to sample.
             :param separate: (True) Separate neutral from non-neutral variants in output.
             :param remove_fixed: (True) Remove sites fixed in the sample from the return value.

             The final sample size is ``2*len(individuals)``.

             The output is a list of tuples.  When separate is False,
             Each tuple is (pair, state), where state is 
             a string encoded as 0 = ancestral and 1 = derived.   Adjacent characters
             in each string are diploid genotypes.  The order of diploids is constant
             across variable sites.

             When ``separate`` is ``True``, the function returns a list of tuples.
             Each tuple is of length two, with the first element being a list of tuples
             for neutral variants following the layout described above.  The second
             element is for selected variants.

			 .. versionadded:: 0.1.4
             )delim");

    py::class_<fwdpy11::singlepop_gm_vec_t,
               singlepop_generalmut_vec_sugar_base>(m,
                                                    "SlocusPopGeneralMutVec")
        .def(py::init<unsigned>(), py::arg("N"),
             "Construct object with N diploids.")
        .def(py::init<const fwdpy11::singlepop_gm_vec_t::dipvector_t&,
                      const fwdpy11::singlepop_gm_vec_t::gcont_t&,
                      const fwdpy11::singlepop_gm_vec_t::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
        .def(py::init<const fwdpy11::singlepop_gm_vec_t&>(),
             R"delim(
                Copy constructor.

                .. versionadded:: 0.1.4
                )delim")
        .def_static(
            "create",
            [](fwdpy11::dipvector_t& diploids, fwdpy11::gcont_t& gametes,
               std::vector<KTfwd::generalmut_vec>& mutations, py::tuple args) {
                if (args.size() == 0)
                    {
                        return fwdpy11::singlepop_gm_vec_t::create(
                            diploids, gametes, mutations);
                    }
                auto& fixations
                    = args[0].cast<std::vector<KTfwd::generalmut_vec>&>();
                auto& fixation_times
                    = args[1].cast<std::vector<KTfwd::uint_t>&>();
                auto g = args[2].cast<KTfwd::uint_t>();
                return fwdpy11::singlepop_gm_vec_t::create_with_fixations(
                    diploids, gametes, mutations, fixations, fixation_times,
                    g);
            })
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
                      ":class:`fwdpy11.VecGeneralMutVec`.")
        .def_readonly("gametes", &fwdpy11::singlepop_gm_vec_t::gametes,
                      GAMETES_DOCSTRING)
        .def_readonly("mcounts", &fwdpy11::singlepop_gm_vec_t::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly("fixations", &fwdpy11::singlepop_gm_vec_t::fixations,
                      "A list of :class:`fwdpy11.VecGeneralMutVec`.")
        .def_readonly("fixation_times",
                      &fwdpy11::singlepop_gm_vec_t::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def_readonly("popdata", &fwdpy11::singlepop_gm_vec_t::popdata,
                      POPDATA_DOCSTRING)
        .def_readwrite("popdata_user",
                       &fwdpy11::singlepop_gm_vec_t::popdata_user,
                       POPDATA_USER_DOCSTRING)
        .def(py::pickle(

            [](const fwdpy11::singlepop_gm_vec_t& pop) -> py::object {
                auto pb = py::bytes(pop.serialize());
                py::list pdata;
                for (auto& d : pop.diploids)
                    {
                        pdata.append(d.parental_data);
                    }
                return py::make_tuple(std::move(pb), std::move(pdata),
                                      pop.popdata, pop.popdata_user);
            },
            [](py::object pickled) {
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 4)
                    {
                        throw std::runtime_error(
                            "expected tuple with 4 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto l = t[1].cast<py::list>();
                auto rv = std::unique_ptr<fwdpy11::singlepop_gm_vec_t>(
                    new fwdpy11::singlepop_gm_vec_t(s));
                for (std::size_t i = 0; i < rv->diploids.size(); ++i)
                    {
                        rv->diploids[i].parental_data = l[i];
                    }
                rv->popdata = t[2];
                rv->popdata_user = t[3];
                return rv;
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
                         return py::make_tuple(s.first, s.second);
                     }
                 py::list rv = py::cast(
                     KTfwd::sample(rng.get(), pop, nsam, remove_fixed));
                 return rv;
             },
             py::arg("rng"), py::arg("nsam"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample random diploids *with replacement*.
             
             :param rng: A :class:`fwdpy11.GSLrng`
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
                :func:`~fwdpy11.SlocusPopGeneralMutVec.sample_ind`.

			 .. versionadded:: 0.1.4
             )delim")
        .def("sample_ind",
             [](const fwdpy11::singlepop_gm_vec_t& pop,
                const std::vector<std::size_t>& individuals,
                const bool separate, const bool remove_fixed) -> py::object {
                 if (separate)
                     {
                         auto s = KTfwd::sample_separate(pop, individuals,
                                                         remove_fixed);
                         return py::make_tuple(s.first, s.second);
                     }
                 py::list rv
                     = py::cast(KTfwd::sample(pop, individuals, remove_fixed));
                 return rv;
             },
             py::arg("individuals"), py::arg("separate") = true,
             py::arg("remove_fixed") = true,
             R"delim(
             Sample specific diploids.
             
             :param rng: A :class:`fwdpy11.GSLrng`
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

			 .. versionadded:: 0.1.4
             )delim");
}
