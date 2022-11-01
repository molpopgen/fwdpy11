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
// This file is copied from the fwdpp version 0.7.4 examples folder.
// Author is Kevin Thornton.
#ifndef FWDPY11_EVOLVE_GENERATION_TS
#define FWDPY11_EVOLVE_GENERATION_TS

#include <cstdint>
#include <algorithm>
#include <memory>
#include <vector>
#include <tuple>
#include <gsl/gsl_randist.h>

#include <fwdpp/util.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/recording/diploid_offspring.hpp>
#include <fwdpp/ts/recording/edge_buffer.hpp>
#include <fwdpp/ts/recording/mutations.hpp>
#include <fwdpy11/discrete_demography/simulation.hpp>

namespace fwdpy11
{

    template <typename poptype, typename rng_t, typename genetic_param_holder>
    std::pair<fwdpp::ts::mut_rec_intermediates, fwdpp::ts::mut_rec_intermediates>
    generate_offspring(const rng_t& rng,
                       const std::pair<std::size_t, std::size_t> parent_indexes,
                       poptype& pop, typename poptype::diploid_t& offspring,
                       genetic_param_holder& genetics)
    {
        auto offspring_data = fwdpp::ts::generate_offspring(
            rng.get(), parent_indexes, fwdpp::ts::selected_variants_only(), pop,
            genetics, offspring);
#ifndef NDEBUG
        for (auto& m : offspring_data.first.mutation_keys)
            {
                auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
                assert(std::distance(itr.first, itr.second) == 1);
            }
        for (auto& m : offspring_data.second.mutation_keys)
            {
                auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
                assert(std::distance(itr.first, itr.second) == 1);
            }
#endif
        return offspring_data;
    }

    inline std::pair<fwdpp::ts::table_index_t, fwdpp::ts::table_index_t>
    parent_nodes_from_metadata(const std::size_t i,
                               const std::vector<fwdpy11::DiploidMetadata>& metadata,
                               const bool swapped)
    {
        auto rv = std::make_pair(metadata[i].nodes[0], metadata[i].nodes[1]);
        if (swapped)
            {
                std::swap(rv.first, rv.second);
            }
        return rv;
    }

    template <typename rng_t, typename poptype, typename genetic_param_holder>
    void
    evolve_generation_ts(
        const rng_t& rng, poptype& pop, genetic_param_holder& genetics,
        const fwdpy11::discrete_demography::DiscreteDemographyState&
            current_demographic_state,
        const fwdpy11::discrete_demography::multideme_fitness_lookups<std::uint32_t>&
            fitness_lookup,
        const fwdpy11::discrete_demography::migration_lookup& miglookup,
        const fwdpp::uint_t generation, fwdpp::ts::edge_buffer& new_edge_buffer,
        std::vector<fwdpy11::DiploidGenotype>& offspring,
        std::vector<fwdpy11::DiploidMetadata>& offspring_metadata,
        std::int32_t next_index)
    {
        fwdpp::debug::all_haploid_genomes_extant(pop);

        genetics.haploid_genome_recycling_bin
            = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);

        fwdpp::zero_out_haploid_genomes(pop);

        offspring.clear();
        offspring_metadata.clear();

        // Generate the offspring
        auto next_index_local = next_index;
        using maxdeme_type = typename std::remove_const<
            decltype(current_demographic_state.maxdemes)>::type;
        for (maxdeme_type deme = 0; deme < current_demographic_state.maxdemes; ++deme)
            {
                auto next_N_deme = current_demographic_state.current_deme_parameters
                                       .next_deme_sizes.get()[deme];
                for (decltype(next_N_deme) ind = 0; ind < next_N_deme; ++ind)
                    {
                        // Get the parents
                        auto pdata = fwdpy11::discrete_demography::pick_parents(
                            rng, deme, miglookup,
                            current_demographic_state.current_deme_parameters
                                .current_deme_sizes,
                            current_demographic_state.current_deme_parameters
                                .selfing_rates,
                            current_demographic_state.fitness_bookmark, fitness_lookup);
                        fwdpy11::DiploidGenotype dip{
                            std::numeric_limits<std::size_t>::max(),
                            std::numeric_limits<std::size_t>::max()};
                        //auto p1 = pick1();
                        //auto p2 = pick2(p1);
                        auto offspring_data = generate_offspring(
                            rng, std::make_pair(pdata.parent1, pdata.parent2), pop, dip,
                            genetics);
                        auto p1id = parent_nodes_from_metadata(
                            pdata.parent1, pop.diploid_metadata,
                            offspring_data.first.swapped);
                        auto p2id = parent_nodes_from_metadata(
                            pdata.parent2, pop.diploid_metadata,
                            offspring_data.second.swapped);
                        fwdpp::ts::table_index_t offspring_node_1
                            = fwdpp::ts::record_diploid_offspring(
                                offspring_data.first.breakpoints, p1id, deme, generation,
                                *pop.tables, new_edge_buffer);
                        fwdpp::ts::record_mutations_infinite_sites(
                            offspring_node_1, pop.mutations,
                            offspring_data.first.mutation_keys, *pop.tables);
                        fwdpp::ts::table_index_t offspring_node_2
                            = fwdpp::ts::record_diploid_offspring(
                                offspring_data.second.breakpoints, p2id, deme,
                                generation, *pop.tables, new_edge_buffer);
                        fwdpp::ts::record_mutations_infinite_sites(
                            offspring_node_2, pop.mutations,
                            offspring_data.second.mutation_keys, *pop.tables);

                        // Add metadata for the offspring
                        offspring_metadata.emplace_back(fwdpy11::DiploidMetadata{
                            0.0,
                            0.0,
                            1.,
                            {0, 0, 0},
                            offspring_metadata.size(),
                            {pdata.parent1, pdata.parent2},
                            deme,
                            0,
                            {offspring_node_1, offspring_node_2}});
                        offspring.emplace_back(std::move(dip));

                        next_index_local = offspring_node_2;
                    }
            }
        assert(static_cast<std::size_t>(next_index_local)
               == pop.tables->num_nodes() - 1);
        if (next_index_local
            != static_cast<decltype(next_index_local)>(pop.tables->num_nodes() - 1))
            {
                throw std::runtime_error("error in book-keeping offspring nodes");
            }
    }
} // namespace fwdpy11
#endif
