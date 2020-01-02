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
#include <sstream>
#include <fstream>
#include <type_traits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/types/create_pops.hpp>
#include <fwdpy11/serialization.hpp>
#include <fwdpy11/serialization/Mutation.hpp>
#include <fwdpy11/serialization/Diploid.hpp>
#include "get_individuals.hpp"

namespace py = pybind11;

namespace
{
    static const auto DIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.DiploidVector`.
   )delim";

} // namespace

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidGenotype>);
PYBIND11_MAKE_OPAQUE(fwdpy11::DiploidPopulation::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::DiploidPopulation::mcont_t);

fwdpy11::DiploidPopulation
create_DiploidPopulation_from_tree_sequence(py::object ts);

void
init_DiploidPopulation(py::module& m)
{
    py::class_<fwdpy11::DiploidPopulation, fwdpy11::Population>(
        m, "DiploidPopulation", "Representation of a diploid population")
        .def(py::init<fwdpp::uint_t, double>(),
             "Construct with an unsigned integer "
             "representing the initial "
             "population size.",
             py::arg("N"),
             py::arg("length") = std::numeric_limits<double>::max())
        .def(py::init<const fwdpy11::DiploidPopulation::dipvector_t&,
                      const fwdpy11::DiploidPopulation::gcont_t&,
                      const fwdpy11::DiploidPopulation::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, haploid_genomes, mutations).
             
             .. versionadded:: 0.1.4
             )delim",
             py::arg("diploids"), py::arg("haploid_genomes"),
             py::arg("mutations"))
        .def(py::init<const fwdpy11::DiploidPopulation&>(),
             R"delim(
                Copy constructor

                .. versionadded:: 0.1.4
                )delim")
        .def("clear", &fwdpy11::DiploidPopulation::clear,
             "Clears all population data.")
        .def("__eq__",
             [](const fwdpy11::DiploidPopulation& lhs,
                const fwdpy11::DiploidPopulation& rhs) { return lhs == rhs; })
        .def_readonly("diploids", &fwdpy11::DiploidPopulation::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_static(
            "create",
            [](fwdpy11::DiploidPopulation::dipvector_t& diploids,
               fwdpy11::DiploidPopulation::gcont_t& haploid_genomes,
               fwdpy11::DiploidPopulation::mcont_t& mutations,
               py::args args) -> fwdpy11::DiploidPopulation {
                if (args.size() == 0)
                    {
                        return fwdpy11::create_wrapper<
                            fwdpy11::DiploidPopulation>()(
                            std::move(diploids), std::move(haploid_genomes),
                            std::move(mutations));
                    }
                auto& fixations
                    = args[0].cast<fwdpy11::DiploidPopulation::mcont_t&>();
                auto ftimes = args[1].cast<std::vector<fwdpp::uint_t>>();
                auto g = args[2].cast<fwdpp::uint_t>();
                return fwdpy11::create_wrapper<fwdpy11::DiploidPopulation>()(
                    std::move(diploids), std::move(haploid_genomes),
                    std::move(mutations), std::move(fixations),
                    std::move(ftimes), g);
            },
            py::arg("diploids"), py::arg("haploid_genomes"),
            py::arg("mutations"),
            R"delim(
        Create a new object from input data.
        Unlike the constructor method, this method results
        in no temporary copies of input data.

        :param diplods: A :class:`fwdpy11.VecDiploid`
        :param haploid_genomes: A :class:`fwdpy11.VecGamete`
        :param mutations: A :class:`fwdpy11.VecMutation`
        :param args: Fixations, fixation times, and generation

        :rtype: :class:`fwdpy11.DiploidPopulation`

        See :ref:`popobjects` for example use.

        When passing in extra args, they must be the following:

        fixations: A :class:`fwdpy11.VecMutation`
        fixation times: A :class:`fwdpy11.VecUint32`
        generation: A non-negative integer

        It is required that len(fixations) == len(fixation times).

        The result of passing in these extra args will be an object
        with its fixation data populated and its generation set
        to the input value.

        .. versionadded:: 0.1.4
        )delim")
        .def(py::pickle(
            [](const fwdpy11::DiploidPopulation& pop) -> py::object {
                std::ostringstream o;
                fwdpy11::serialization::serialize_details(o, &pop);
                auto pb = py::bytes(o.str());
                return py::object(std::move(pb));
            },
            [](py::object pickled) -> fwdpy11::DiploidPopulation {
                auto s = pickled.cast<py::bytes>().cast<std::string>();
                fwdpy11::DiploidPopulation pop(
                    1, std::numeric_limits<double>::max());
                std::istringstream in(std::move(s));
                fwdpy11::serialization::deserialize_details()(in, pop);
                return pop;
            }))
        .def(
            "sample",
            [](const fwdpy11::DiploidPopulation& pop,
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
        .def(
            "sample",
            [](const fwdpy11::DiploidPopulation& pop,
               const fwdpy11::GSLrng_t& rng, const std::uint32_t nsam,
               const bool haplotype, const bool remove_fixed) {
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
        .def("add_mutations", &fwdpy11::DiploidPopulation::add_mutations)
        .def(
            "dump_to_file",
            [](const fwdpy11::DiploidPopulation& pop,
               const std::string filename) {
                std::ofstream out(filename.c_str(), std::ios_base::binary);
                if (!out)
                    {
                        throw std::runtime_error(
                            "could not open file for writing");
                    }
                fwdpy11::serialization::serialize_details(out, &pop);
                out.close();
            },
            "Write a population to a file in binary format.")
        .def_static(
            "load_from_file",
            [](const std::string filename) {
                std::ifstream in(filename.c_str(), std::ios_base::binary);
                if (!in)
                    {
                        throw std::runtime_error(
                            "could not open file for reading");
                    }
                fwdpy11::DiploidPopulation pop(
                    1, std::numeric_limits<double>::max());
                fwdpy11::serialization::deserialize_details()(in, pop);
                return pop;
            },
            "Load a population from a binary file.")
        .def(
            "pickle_to_file",
            [](const fwdpy11::DiploidPopulation& self, py::object f) {
                auto dump = py::module::import("pickle").attr("dump");
                dump(py::make_tuple(
                         self.diploids.size(), self.haploid_genomes.size(),
                         self.mutations.size(), self.fixations.size(),
                         self.generation, self.tables.genome_length()),
                     f);
                for (auto& d : self.diploids)
                    {
                        dump(d, f);
                    }
                for (auto& g : self.haploid_genomes)
                    {
                        dump(g, f);
                    }
                for (auto& m : self.mutations)
                    {
                        dump(m, f);
                    }
                for (auto& m : self.fixations)
                    {
                        dump(m, f);
                    }
                dump(self.fixation_times, f);
                dump(self.mcounts, f);
                dump(self.mcounts_from_preserved_nodes, f);
                dump(py::make_tuple(self.diploid_metadata.size(),
                                    self.ancient_sample_metadata.size()),
                     f);
                for (auto& md : self.diploid_metadata)
                    {
                        dump(md, f);
                    }
                for (auto& md : self.ancient_sample_metadata)
                    {
                        dump(md, f);
                    }
                dump(py::make_tuple(self.tables.node_table.size(),
                                    self.tables.edge_table.size(),
                                    self.tables.mutation_table.size(),
                                    self.tables.site_table.size()),
                     f);
                for (auto& n : self.tables.node_table)
                    {
                        dump(n, f);
                    }
                for (auto& e : self.tables.edge_table)
                    {
                        dump(e, f);
                    }
                for (auto& m : self.tables.mutation_table)
                    {
                        dump(m, f);
                    }
                for (auto& s : self.tables.site_table)
                    {
                        dump(s, f);
                    }
                dump(self.tables.preserved_nodes, f);
            },
            R"delim(
             Pickle the population to an open file.

             This function may be preferred over 
             the direct pickling method because it uses less
             memory.  It is, however, slower.

             To read the population back in, you must call
             :func:`fwdpy11.DiploidPopulation.load_from_pickle_file`.

             :param f: A handle to an open file

             .. versionadded:: 0.3.0
             )delim")
        .def_static(
            "load_from_pickle_file",
            [](py::object f) {
                auto load = py::module::import("pickle").attr("load");
                py::tuple popdata = load(f);
                fwdpy11::DiploidPopulation rv(popdata[0].cast<fwdpp::uint_t>(),
                                              popdata[5].cast<double>());
                rv.generation
                    = popdata[4]
                          .cast<decltype(
                              fwdpy11::DiploidPopulation::generation)>();
                auto ndips = popdata[0].cast<std::size_t>();
                auto ngams = popdata[1].cast<std::size_t>();
                auto nmuts = popdata[2].cast<std::size_t>();
                auto nfixations = popdata[3].cast<std::size_t>();
                rv.diploids.clear();
                rv.haploid_genomes.clear();
                rv.mutations.clear();
                rv.fixations.clear();
                rv.diploids.reserve(ndips);
                for (std::size_t i = 0; i < ndips; ++i)
                    {
                        rv.diploids.push_back(
                            load(f).cast<fwdpy11::DiploidGenotype>());
                    }
                rv.haploid_genomes.reserve(ngams);
                for (std::size_t i = 0; i < ngams; ++i)
                    {
                        rv.haploid_genomes.push_back(
                            load(f).cast<fwdpp::haploid_genome>());
                    }
                rv.mutations.reserve(nmuts);
                for (std::size_t i = 0; i < nmuts; ++i)
                    {
                        rv.mutations.push_back(
                            load(f).cast<fwdpy11::Mutation>());
                    }
                rv.fixations.reserve(nfixations);
                for (std::size_t i = 0; i < nfixations; ++i)
                    {
                        rv.fixations.push_back(
                            load(f).cast<fwdpy11::Mutation>());
                    }
                rv.fixation_times
                    = load(f).cast<decltype(rv.fixation_times)>();
                rv.mcounts = load(f).cast<decltype(rv.mcounts)>();
                rv.mcounts_from_preserved_nodes
                    = load(f)
                          .cast<decltype(rv.mcounts_from_preserved_nodes)>();
                py::tuple metadata_data = load(f);
                rv.diploid_metadata.clear();
                rv.ancient_sample_metadata.clear();
                auto lmd = metadata_data[0].cast<std::size_t>();
                auto lamd = metadata_data[1].cast<std::size_t>();
                rv.diploid_metadata.reserve(lmd);
                for (std::size_t i = 0; i < lmd; ++i)
                    {
                        rv.diploid_metadata.push_back(
                            load(f).cast<fwdpy11::DiploidMetadata>());
                    }
                rv.ancient_sample_metadata.reserve(lamd);
                for (std::size_t i = 0; i < lamd; ++i)
                    {
                        rv.ancient_sample_metadata.push_back(
                            load(f).cast<fwdpy11::DiploidMetadata>());
                    }
                py::tuple table_data = load(f);
                auto table_len = table_data[0].cast<std::size_t>();
                rv.tables.clear();
                rv.tables.node_table.reserve(table_len);
                for (std::size_t i = 0; i < table_len; ++i)
                    {
                        rv.tables.node_table.push_back(
                            load(f).cast<fwdpp::ts::node>());
                    }
                table_len = table_data[1].cast<std::size_t>();
                rv.tables.edge_table.reserve(table_len);
                for (std::size_t i = 0; i < table_len; ++i)
                    {
                        rv.tables.edge_table.push_back(
                            load(f).cast<fwdpp::ts::edge>());
                    }
                table_len = table_data[2].cast<std::size_t>();
                rv.tables.mutation_table.reserve(table_len);
                for (std::size_t i = 0; i < table_len; ++i)
                    {
                        rv.tables.mutation_table.push_back(
                            load(f).cast<fwdpp::ts::mutation_record>());
                    }
                table_len = table_data[3].cast<std::size_t>();
                rv.tables.site_table.reserve(table_len);
                for (std::size_t i = 0; i < table_len; ++i)
                    {
                        rv.tables.site_table.push_back(
                            load(f).cast<fwdpp::ts::site>());
                    }

                rv.tables.preserved_nodes
                    = load(f).cast<decltype(rv.tables.preserved_nodes)>();
                rv.tables.build_indexes();
                return rv;
            },
            R"delim(
            Read in a pickled population from a file.
            The file muse have been generated by
            a call to :func:`fwdpy11.DiploidPopulation.pickle_to_file`.

            :param f: A handle to a file opened in 'rb' mode.

            .. versionadded: 0.3.0
            )delim")
        .def_static(
            "create_from_tskit",
            [](py::object ts) {
                return create_DiploidPopulation_from_tree_sequence(ts);
            },
            R"delim(
        Create a new object from an tskit.TreeSequence

        :param ts: A tree sequence from tskit
        :type ts: tskit.TreeSequence

        :rtype: :class:`fwdpy11.DiploidPopulation`
        :returns: A population object with an initialized
        :class:`fwdpy11.TableCollection`

        .. versionadded:: 0.2.0

        .. note::

            In general, initializing a population using
            the output from a coalescent simulation is
            a tricky business.  There are issues of
            parameter scaling and the appropriateness
            of the coalescent model itself. A key issue
            is that your input tree sequence must have
            node times in the correct time units! (Generations,
            for example.) See :ref:`ts` for more discussion
        )delim");
}
