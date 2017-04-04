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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types.hpp>

namespace py = pybind11;

using fwdpp_popgenmut_base = fwdpy11::singlepop_t::popbase_t;
using singlepop_sugar_base = fwdpy11::singlepop_t::base;
using multilocus_sugar_base = fwdpy11::multilocus_t::base;
using multilocus_popgenmut_base = multilocus_sugar_base::popbase_t;
using singlepop_generalmut_vec_sugar_base = fwdpy11::singlepop_gm_vec_t::base;
using singlepop_generalmut_vec_base
    = singlepop_generalmut_vec_sugar_base::popbase_t;

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

PYBIND11_PLUGIN(fwdpy11_types)
{
    py::module m("fwdpy11_types", "Wrap C++ types specific to fwdpy11.");

    py::class_<fwdpy11::GSLrng_t>(
        m, "GSLrng", "Random number generator based on a mersenne twister.")
        .def(py::init<unsigned>(),
             "Constructor takes unsigned integer as a seed");

    py::class_<fwdpy11::diploid_t>(
        m, "SingleLocusDiploid",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_readonly("first", &fwdpy11::diploid_t::first,
                      "Key to first gamete. (read-only)")
        .def_readonly("second", &fwdpy11::diploid_t::second,
                      "Key to second gamete. (read-onle)")
        .def_readonly("w", &fwdpy11::diploid_t::w, "Fitness. (read-only)")
        .def_readonly("g", &fwdpy11::diploid_t::g,
                      "Genetic value (read-only).")
        .def_readonly("e", &fwdpy11::diploid_t::e,
                      "Random/environmental effects (read-only).")
        .def("__getstate__",
             [](const fwdpy11::diploid_t& d) {
                 return py::make_tuple(d.first, d.second, d.w, d.g, d.e);
             })
        .def("__setstate__", [](fwdpy11::diploid_t& d, py::tuple t) {
            new (&d) fwdpy11::diploid_t(t[0].cast<std::size_t>(),
                                        t[1].cast<std::size_t>());
            d.w = t[2].cast<double>();
            d.g = t[3].cast<double>();
            d.e = t[4].cast<double>();
        });

    py::bind_vector<fwdpy11::dipvector_t>(
        m, "DiploidContainer",
        "C++ representation of a list of "
        ":class:`fwdpy11.fwdpy11_types."
        "SingleLocusDiploid`.  Typically, access will be read-only.");
    py::bind_vector<std::vector<KTfwd::uint_t>>(m, "VectorUint32");
    py::bind_vector<fwdpy11::gcont_t>(m, "GameteContainer",
                                      "C++ representations of a list of "
                                      ":class:`fwdpy11.fwdpp_types.Gamete`.  "
                                      "Typically, access will be read-only.");
    py::bind_vector<fwdpy11::mcont_t>(
        m, "MutationContainer", "C++ representation of a list of "
                                ":class:`fwdpy11.fwdpp_types.Mutation`.  "
                                "Typically, access will be read-only.");

    // expose the base classes for population types
    py::class_<fwdpp_popgenmut_base>(m, "SpopMutationBase");
    py::class_<multilocus_popgenmut_base>(m, "MlocusMutationBase");
    py::class_<singlepop_sugar_base, fwdpp_popgenmut_base>(m, "SinglepopBase");
    py::class_<multilocus_sugar_base, multilocus_popgenmut_base>(m,
                                                                 "MlocusBase");

    py::class_<singlepop_generalmut_vec_base>(m, "SpopGeneralMutVecBase");
    py::class_<singlepop_generalmut_vec_sugar_base,
               singlepop_generalmut_vec_base>(m, "SpopGeneralMutVecSugarBase");

    // Expose the type based on fwdpp's "sugar" layer
    py::class_<fwdpy11::singlepop_t, singlepop_sugar_base>(
        m, "Spop", "Population object representing a single deme.")
        .def(py::init<unsigned>(),
             "Construct with an unsigned integer representing the initial "
             "population size.")
        .def("clear", &fwdpy11::singlepop_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::singlepop_t::generation,
                      R"delim(
                      The current generation. A population starts at 
                      generation 0:

                        .. testcode:: 

                            import fwdpy11
                            p = fwdpy11.Spop(1000)
                            print(p.generation)
                            import fwdpy11.wright_fisher as wf
                            p = wf.quick_sim(100)
                            print(p.generation)

                        The output is:

                        .. testoutput::

                          0
                          100
                        )delim")
        .def_readonly("N", &fwdpy11::singlepop_t::N,
                      R"delim(
                      The current population size.

                      .. testcode:: 

                          import fwdpy11
                          p = fwdpy11.Spop(1000)
                          print(p.N)
                          import numpy as np
                          import fwdpy11.wright_fisher as wf
                          #Evolve to a final N of 500
                          nlist = np.array([p.N]*100 + [p.N/2]*100,dtype = np.uint32)
                          rng=fwdpy11.GSLrng(101)
                          wf.evolve(rng,p,nlist)
                          print(p.N)
                          print(p.generation)

                      The output from the above is:

                      .. testoutput::
                      
                        1000
                        500
                        200
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
        .def("__getstate__",
             [](const fwdpy11::singlepop_t& pop) {
                 return py::bytes(pop.serialize());
             })
        .def("__setstate__",
             [](fwdpy11::singlepop_t& p, py::bytes s) {
                 new (&p) fwdpy11::singlepop_t(s);
             })
        .def("__eq__",
             [](const fwdpy11::singlepop_t& lhs,
                const fwdpy11::singlepop_t& rhs) { return lhs == rhs; });

    py::class_<fwdpy11::multilocus_t, multilocus_sugar_base>(m, "MlocusPop")
        .def(py::init<unsigned, unsigned>(), py::arg("N"), py::arg("nloci"),
                "Construct with population size and number of loci.")
        .def("clear", &fwdpy11::multilocus_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::multilocus_t::generation,
                      "The current generation.")
        .def_readonly("N", &fwdpy11::multilocus_t::N,
                      "Curent population size.")
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
        .def("__getstate__",
             [](const fwdpy11::multilocus_t& pop) {
                 return py::bytes(pop.serialize());
             })
        .def("__setstate__",
             [](fwdpy11::multilocus_t& p, py::bytes s) {
                 new (&p) fwdpy11::multilocus_t(s);
             })
        .def("__eq__",
             [](const fwdpy11::multilocus_t& lhs,
                const fwdpy11::multilocus_t& rhs) { return lhs == rhs; });

    py::class_<fwdpy11::singlepop_gm_vec_t,
               singlepop_generalmut_vec_sugar_base>(
        m, "SpopGeneralMutVec", "Single-deme object using "
                                ":class:`fwpy11.fwdpp_types.GeneralMutVec` as "
                                "the mutation type.")
        .def(py::init<unsigned>(), py::arg("N"),
                "Construct object with N diploids.")
        .def("clear", &fwdpy11::singlepop_gm_vec_t::clear,
             "Clears all population data.")
        .def_readonly("generation", &fwdpy11::singlepop_gm_vec_t::generation,
                      "The current generation.")
        .def_readonly("N", &fwdpy11::singlepop_gm_vec_t::N,
                      "Curent population size.")
        .def_readonly("diploids", &fwdpy11::singlepop_gm_vec_t::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_readonly(
            "mutations", &fwdpy11::singlepop_gm_vec_t::mutations,
            "A list of :class:`fwdpy11.fwdpp_types.VectorGeneralMutVec`.")
        .def_readonly("gametes", &fwdpy11::singlepop_gm_vec_t::gametes,
                      GAMETES_DOCSTRING)
        .def_readonly("mcounts", &fwdpy11::singlepop_gm_vec_t::mcounts,
                      MCOUNTS_DOCSTRING)
        .def_readonly(
            "fixations", &fwdpy11::singlepop_gm_vec_t::fixations,
            "A list of :class:`fwdpy11.fwdpp_types.VectorGeneralMutVec`.")
        .def_readonly("fixation_times",
                      &fwdpy11::singlepop_gm_vec_t::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def("__getstate__",
             [](const fwdpy11::singlepop_gm_vec_t& pop) {
                 return py::bytes(pop.serialize());
             })
        .def("__setstate__",
             [](fwdpy11::singlepop_gm_vec_t& p, py::bytes s) {
                 new (&p) fwdpy11::singlepop_gm_vec_t(s);
             })
        .def("__eq__", [](const fwdpy11::singlepop_gm_vec_t& lhs,
                          const fwdpy11::singlepop_gm_vec_t& rhs) {
            return lhs == rhs;
        });
    return m.ptr();
}
