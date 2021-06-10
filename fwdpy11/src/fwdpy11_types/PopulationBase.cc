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
    std::int32_t g;
    decltype(fwdpy11::Mutation::xtra) label;
    std::int16_t neutral;
};

namespace
{
    py::array
    make_flattened_Mutation_array(
        const fwdpy11::Population::mutation_container& mutations)
    {
        std::vector<flattened_Mutation> vfm;
        vfm.reserve(mutations.size());
        for (auto&& m : mutations)
            {
                vfm.push_back(
                    flattened_Mutation{m.pos, m.s, m.h, m.g, m.xtra, m.neutral});
            }
        return fwdpy11::make_1d_array_with_capsule(std::move(vfm));
    }
} // namespace

PYBIND11_MAKE_OPAQUE(fwdpy11::Population::genome_container);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::mutation_container);

void
init_PopulationBase(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(flattened_Mutation, pos, s, h, g, label, neutral);

    py::class_<fwdpy11::Population>(m, "PopulationBase")
        .def_readonly("_N", &fwdpy11::Population::N)
        .def_readonly("_generation", &fwdpy11::Population::generation)
        .def_readonly("_is_simulating", &fwdpy11::Population::is_simulating)
        .def_readonly("_mutations", &fwdpy11::Population::mutations)
        .def_property_readonly("_mutations_ndarray",
                               [](const fwdpy11::Population& self) {
                                   return make_flattened_Mutation_array(self.mutations);
                               })
        .def_property_readonly("_mcounts",
                               [](const fwdpy11::Population& self) {
                                   return fwdpy11::make_1d_ndarray_readonly(
                                       self.mcounts);
                               })
        .def_property_readonly("_mcounts_ancient_samples",
                               [](const fwdpy11::Population& self) {
                                   return fwdpy11::make_1d_ndarray_readonly(
                                       self.mcounts_from_preserved_nodes);
                               })
        .def_property_readonly(
            "_mut_lookup",
            [](const fwdpy11::Population& pop) -> py::object {
                py::dict rv;
                for (std::size_t i = 0; i < pop.mutations.size(); ++i)
                    {
                        if (pop.mcounts[i] > 0
                            || (!pop.mcounts_from_preserved_nodes.empty()
                                && pop.mcounts_from_preserved_nodes[i] > 0))
                            {
                                auto pos_handle = py::cast(pop.mutations[i].pos);
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
            })
        .def(
            "_mutation_indexes",
            [](const fwdpy11::Population& pop, const double pos) -> py::object {
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
            py::arg("pos"))
        .def_readonly("_haploid_genomes", &fwdpy11::Population::haploid_genomes)
        .def_readonly("_fixations", &fwdpy11::Population::fixations)
        .def_readonly("_fixation_times", &fwdpy11::Population::fixation_times)
        .def("_find_mutation_by_key",
             [](const fwdpy11::Population& pop,
                const std::tuple<double, double, fwdpp::uint_t>& key,
                const std::int64_t offset) {
                 return pop.find_mutation_by_key(key, offset);
             })
        .def("_find_fixation_by_key",
             [](const fwdpy11::Population& pop,
                const std::tuple<double, double, fwdpp::uint_t>& key,
                const std::int64_t offset) {
                 return pop.find_fixation_by_key(key, offset);
             })
        // TODO: why does readwrite fail?
        .def_readonly("_tables", &fwdpy11::Population::tables)
        .def_property_readonly("_genetic_values",
                               [](const fwdpy11::Population& self) {
                                   return fwdpy11::make_2d_ndarray_readonly(
                                       self.genetic_value_matrix, self.N,
                                       self.genetic_value_matrix.size() / self.N);
                               })
        .def_property_readonly(
            "_ancient_sample_genetic_values", [](const fwdpy11::Population& self) {
                if (self.ancient_sample_metadata_size() == 0
                    || self.ancient_sample_genetic_value_matrix.empty())
                    {
                        return fwdpy11::make_1d_ndarray_readonly(
                            self.ancient_sample_genetic_value_matrix);
                    }
                return fwdpy11::make_2d_ndarray_readonly(
                    self.ancient_sample_genetic_value_matrix,
                    self.ancient_sample_metadata_size(),
                    self.ancient_sample_genetic_value_matrix.size()
                        / self.ancient_sample_metadata_size());
            });
}
