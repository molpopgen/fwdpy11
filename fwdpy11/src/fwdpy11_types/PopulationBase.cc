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
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

struct flattened_Mutation
{
    double pos, s, h;
    fwdpp::uint_t g;
    decltype(fwdpy11::Mutation::xtra) label;
    std::int16_t neutral;
};

namespace
{
    static const auto MCOUNTS_DOCSTRING = R"delim(
    List of number of occurrences of elements in 
    a population objecst "mutations" container.

    The values are unsigned 32-bit integers.  

    .. note::
        Some values may be 0.  These represent *extinct* variants.  You will typically want to avoid processing such mutations.
)delim";

    static const auto FIXATIONS_DOCSTRING
        = R"delim(A :class:`fwdpy11.VecMutation` of fixed variants.)delim";

    static const auto FIXATION_TIMES_DOCSTRING
        = R"delim(A list of fixation times corresponding to the elements in "fixations" for this type.)delim";

    static const auto GAMETES_DOCSTRING
        = R"delim(A :class:`fwdpy11.HaploidGenomeVector`.)delim";

    static const auto MUTATIONS_DOCSTRING = R"delim(
    List of :class:`fwdpy11.Mutation`.

    .. note:: 
        This list contains **both** extinct *and* extant mutations.  
        To distinguish them, use the locations of nonzero values in "mcounts" 
        for an instance of this type."
    )delim";

    py::array
    make_flattened_Mutation_array(
        const fwdpy11::Population::mcont_t& mutations)
    {
        std::vector<flattened_Mutation> vfm;
        vfm.reserve(mutations.size());
        for (auto&& m : mutations)
            {
                vfm.push_back(flattened_Mutation{ m.pos, m.s, m.h, m.g, m.xtra,
                                                  m.neutral });
            }
        return fwdpy11::make_1d_array_with_capsule(std::move(vfm));
    }
} // namespace

PYBIND11_MAKE_OPAQUE(fwdpy11::Population::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::mcont_t);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidMetadata>);

void
init_PopulationBase(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(flattened_Mutation, pos, s, h, g, label, neutral);

    py::class_<fwdpy11::Population>(m, "PopulationBase",
                                    "Abstract base class for populations "
                                    "based on :class:`fwdpy11.Mutation`")
        .def_readonly("N", &fwdpy11::Population::N, "Current population size.")
        .def_readonly("generation", &fwdpy11::Population::generation,
                      "Curent generation.")
        .def_readonly("mutations", &fwdpy11::Population::mutations,
                      MUTATIONS_DOCSTRING)
        .def_property_readonly(
            "mutations_ndarray",
            [](const fwdpy11::Population& self) {
                return make_flattened_Mutation_array(self.mutations);
            },
            R"delim(
                               Return readonly numpy.ndarray of mutation data.
                               The data are returned as a copy.

                               The esizes and heffects fields are not part of the 
                               array, as their size is not constant and therefore
                               they canot be part of a numpy "dtype".

                               .. versionadded:: 0.4.0
                               )delim")
        .def_property_readonly(
            "mcounts",
            [](const fwdpy11::Population& self) {
                return fwdpy11::make_1d_ndarray_readonly(self.mcounts);
            },
            MCOUNTS_DOCSTRING)
        .def_property_readonly(
            "mcounts_ancient_samples",
            [](const fwdpy11::Population& self) {
                return fwdpy11::make_1d_ndarray_readonly(
                    self.mcounts_from_preserved_nodes);
            },
            "The contribution that ancient samples make to mutation counts")
        .def_readwrite("diploid_metadata",
                       &fwdpy11::Population::diploid_metadata,
                       "Container of diploid metadata.")
        .def_readwrite("ancient_sample_metadata",
                       &fwdpy11::Population::ancient_sample_metadata,
                       "Container of metadata for ancient samples.")
        .def_property_readonly(
            "mut_lookup",
            [](const fwdpy11::Population& pop) -> py::object {
                py::dict rv;
                for (std::size_t i = 0; i < pop.mutations.size(); ++i)
                    {
                        if (pop.mcounts[i])
                            {
                                auto pos_handle
                                    = py::cast(pop.mutations[i].pos);
                                if (rv.contains(pos_handle))
                                    {
                                        py::list l = rv[pos_handle];
                                        l.append(i);
                                    }
                                else
                                    {
                                        py::list l;
                                        l.append(i);
                                        rv[pos_handle] = l;
                                    }
                            }
                    }
                if (rv.size() == 0)
                    {
                        return py::none();
                    }
                return py::object(std::move(rv));
            },
            R"delim(
            Mutation position lookup table.
            The format is a dict whose keys 
            are mutation positions and values
            are lists of indexes.  If no
            extant mutations are present in the 
            population, the value of this property is
            None.
            )delim")
        .def(
            "mutation_indexes",
            [](const fwdpy11::Population& pop,
               const double pos) -> py::object {
                auto r = pop.mut_lookup.equal_range(pos);
                if (r.first == r.second)
                    {
                        return py::none();
                    }
                std::vector<std::size_t> rv;
                for (auto i = r.first; i != r.second; ++i)
                    {
                        rv.push_back(i->second);
                    }
                std::sort(begin(rv), end(rv));
                return fwdpy11::make_1d_array_with_capsule(std::move(rv));
            },
            R"delim(
             Get indexes associated with a mutation position.
             
             :param pos: A position
             :type pos: float
             :return: Indexes in mutation/mutation counts container associated with pos.
             :rtype: numpy.array
             
             Returns None if pos does not refer to an extant variant.  Otherwise, 
             returns a list.
             )delim",
            py::arg("pos"))
        .def_readonly("haploid_genomes", &fwdpy11::Population::haploid_genomes,
                      GAMETES_DOCSTRING)
        .def_readonly("fixations", &fwdpy11::Population::fixations,
                      FIXATIONS_DOCSTRING)
        .def_readonly("fixation_times", &fwdpy11::Population::fixation_times,
                      FIXATION_TIMES_DOCSTRING)
        .def(
            "find_mutation_by_key",
            [](const fwdpy11::Population& pop,
               const std::tuple<double, double, fwdpp::uint_t>& key,
               const std::int64_t offset) {
                return pop.find_mutation_by_key(key, offset);
            },
            py::arg("pop"), py::arg("offset") = 0,
            R"delim(
             Find a mutation by key.
             
             :param key: A mutation key. See :func:`fwdpy11.Mutation.key`.
             :type key: tuple
             :param offset: Offset to start search in mutation container.
             :type offset: int

             :rtype: int

             :returns: Index of mutation if found, or -1 otherwise.

             .. versionadded:: 0.2.0
             )delim")
        .def(
            "find_fixation_by_key",
            [](const fwdpy11::Population& pop,
               const std::tuple<double, double, fwdpp::uint_t>& key,
               const std::int64_t offset) {
                return pop.find_fixation_by_key(key, offset);
            },
            py::arg("pop"), py::arg("offset") = 0,
            R"delim(
             Find a fixation by key.
             
             :param key: A mutation key. See :func:`fwdpy11.Mutation.key`.
             :type key: tuple
             :param offset: Offset to start search in fixation container.
             :type offset: int

             :rtype: int

             :returns: Index of fixation if found, or -1 otherwise.

             .. versionadded:: 0.2.0
             )delim")
        // TODO: why does readwrite fail?
        .def_readonly("tables", &fwdpy11::Population::tables,
                      R"delim(
                Give access to the population's 
                :class:`fwdpy11.TableCollection`
                )delim")
        .def_property_readonly(
            "genetic_values",
            [](const fwdpy11::Population& self) {
                return fwdpy11::make_2d_ndarray_readonly(
                    self.genetic_value_matrix, self.N,
                    self.genetic_value_matrix.size() / self.N);
            },
            R"delim(
        Return the genetic values as a readonly 2d numpy.ndarray.
        
        Rows are individuals.  Columns are genetic values.
        
        ..  versionadded 0.3.0
        )delim")
        .def_property_readonly(
            "ancient_sample_genetic_values",
            [](const fwdpy11::Population& self) {
                return fwdpy11::make_2d_ndarray_readonly(
                    self.ancient_sample_genetic_value_matrix,
                    self.ancient_sample_metadata.size(),
                    self.ancient_sample_genetic_value_matrix.size()
                        / self.ancient_sample_metadata.size());
            },
            R"delim(
        Return the genetic values for ancient samples as a readonly 2d 
        numpy.ndarray.
        
        Rows are individuals.  Columns are genetic values.
        
        ..  versionadded 0.3.0
        )delim");
}
