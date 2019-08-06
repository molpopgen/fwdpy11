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
#include <vector>
#include <tuple>
#include <gsl/gsl_randist.h>

#include <fwdpp/util.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>

namespace fwdpy11
{

    template <typename poptype, typename rng_t, typename genetic_param_holder>
    std::pair<fwdpp::ts::mut_rec_intermediates,
              fwdpp::ts::mut_rec_intermediates>
    generate_offspring(
        const rng_t& rng,
        const std::pair<std::size_t, std::size_t> parent_indexes, poptype& pop,
        typename poptype::diploid_t& offspring, genetic_param_holder& genetics)
    {
        auto offspring_data = fwdpp::ts::generate_offspring(
            rng.get(), parent_indexes, fwdpp::ts::selected_variants_only(),
            pop, genetics, offspring);
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

    // NOTE: much of what follows is simply copy-pasted from 
    // the inner workings of fwdpp::ts::table_collection.
    inline void
    split_breakpoints_add_edges(
        const std::vector<double>& breakpoints, std::size_t parent,
        const std::tuple<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>&
            parents,
        const fwdpp::ts::TS_NODE_INT next_index, double genome_length,
        std::vector<std::array<fwdpp::ts::edge_vector, 2>>& temp_edges)
    {
        auto p1 = std::get<0>(parents), p2 = std::get<1>(parents);
        std::size_t p1_idx = (p1 < p2) ? 0 : 1;
        std::size_t p2_idx = !p1_idx;
        if (breakpoints.front() != 0.0)
            {
                temp_edges[parent][p1_idx].push_back(fwdpp::ts::edge{
                    0., breakpoints.front(), p1, next_index });
            }
        // TODO: replace with exception via a debug mode
        assert(std::count(begin(breakpoints), end(breakpoints),
                          std::numeric_limits<double>::max())
               == 1);
        assert(breakpoints.back() == std::numeric_limits<double>::max());
        for (unsigned j = 1; j < breakpoints.size(); ++j)
            {
                double a = breakpoints[j - 1];
                double b = (j < breakpoints.size() - 1) ? breakpoints[j]
                                                        : genome_length;
                if (b <= a)
                    {
                        throw std::invalid_argument("right must be > left");
                    }
                if (j % 2 == 0.)
                    {
                        temp_edges[parent][p1_idx].push_back(
                            fwdpp::ts::edge{ a, b, p1, next_index });
                    }
                else
                    {
                        temp_edges[parent][p2_idx].push_back(
                            fwdpp::ts::edge{ a, b, p2, next_index });
                    }
            }
    }

    inline void
    populate_temp_edges(
        const std::vector<double>& breakpoints, std::size_t parent,
        const std::tuple<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>&
            parents,
        fwdpp::ts::TS_NODE_INT next_index, double genome_length,
        std::vector<std::array<fwdpp::ts::edge_vector, 2>>& temp_edges)
    {
        if (breakpoints.empty())
            {
                auto p = std::get<0>(parents);
                std::size_t p_idx = (p < std::get<1>(parents)) ? 0 : 1;
                temp_edges[parent][p_idx].push_back(
                    fwdpp::ts::edge{ 0., genome_length, p, next_index });
                return;
            }
        auto itr = std::adjacent_find(std::begin(breakpoints),
                                      std::end(breakpoints));
        if (itr == std::end(breakpoints))
            {
                split_breakpoints_add_edges(breakpoints, parent, parents,
                                            next_index, genome_length,
                                            temp_edges);
            }
        else
            {
                // Here, we need to reduce the input
                // breakpoints to only those seen
                // an odd number of times.
                // Even numbers of the same breakpoint
                // are "double x-overs" and thus
                // cannot affect the genealogy.
                std::vector<double> odd_breakpoints;
                auto start = breakpoints.begin();
                while (itr < breakpoints.end())
                    {
                        auto not_equal = std::find_if(
                            itr, breakpoints.end(),
                            [itr](const double d) { return d != *itr; });
                        int even = (std::distance(itr, not_equal) % 2 == 0.0);
                        odd_breakpoints.insert(odd_breakpoints.end(), start,
                                               itr + 1 - even);
                        start = not_equal;
                        itr = std::adjacent_find(start, std::end(breakpoints));
                    }
                odd_breakpoints.insert(odd_breakpoints.end(), start,
                                       breakpoints.end());
                split_breakpoints_add_edges(odd_breakpoints, parent, parents,
                                            next_index, genome_length,
                                            temp_edges);
            }
    }

    inline fwdpp::ts::TS_NODE_INT
    register_diploid_offspring(
        const std::vector<double>& breakpoints, std::size_t parent,
        const std::tuple<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>&
            parents,
        const std::int32_t population, const double time,
        fwdpp::ts::table_collection& tables,
        std::vector<std::array<fwdpp::ts::edge_vector, 2>>& temp_edges)
    {
        auto next_index = tables.emplace_back_node(population, time);
        if (next_index >= std::numeric_limits<fwdpp::ts::TS_NODE_INT>::max())
            {
                throw std::invalid_argument("node index too large");
            }
        populate_temp_edges(breakpoints, parent, parents, next_index,
                            tables.genome_length(), temp_edges);
        return next_index;
    }

    template <typename rng_t, typename poptype, typename pick_parent1_fxn,
              typename pick_parent2_fxn, typename offspring_metadata_fxn,
              typename genetic_param_holder>
    void
    evolve_generation_ts(
        const rng_t& rng, poptype& pop,
        std::vector<std::array<fwdpp::ts::edge_vector, 2>>& temp_edges,
        std::vector<std::pair<std::size_t, std::size_t>> & offsets,
        fwdpp::ts::edge_vector& temp_edges2, genetic_param_holder& genetics,
        const fwdpp::uint_t N_next, const pick_parent1_fxn& pick1,
        const pick_parent2_fxn& pick2,
        const offspring_metadata_fxn& update_offspring,
        const fwdpp::uint_t generation, fwdpp::ts::table_collection& tables,
        std::int32_t first_parental_index, std::int32_t next_index)
    {
        fwdpp::debug::all_haploid_genomes_extant(pop);

        genetics.haploid_genome_recycling_bin
            = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);

        fwdpp::zero_out_haploid_genomes(pop);

        decltype(pop.diploids) offspring(N_next);
        decltype(pop.diploid_metadata) offspring_metadata(N_next);

        // Generate the offspring
        auto next_index_local = next_index;
        for (std::size_t next_offspring = 0; next_offspring < offspring.size();
             ++next_offspring)
            {
                auto p1 = pick1();
                auto p2 = pick2(p1);
                auto& dip = offspring[next_offspring];
                auto offspring_data = generate_offspring(
                    rng, std::make_pair(p1, p2), pop, dip, genetics);
                auto p1id = fwdpp::ts::get_parent_ids(
                    first_parental_index, p1, offspring_data.first.swapped);
                auto p2id = fwdpp::ts::get_parent_ids(
                    first_parental_index, p2, offspring_data.second.swapped);
                next_index_local = register_diploid_offspring(
                    offspring_data.first.breakpoints, p1, p1id, 0, generation,
                    pop.tables, temp_edges);
                fwdpp::ts::record_mutations_infinite_sites(
                    next_index_local, pop.mutations,
                    offspring_data.first.mutation_keys, tables);
                next_index_local = register_diploid_offspring(
                    offspring_data.second.breakpoints, p2, p2id, 0, generation,
                    pop.tables, temp_edges);
                fwdpp::ts::record_mutations_infinite_sites(
                    next_index_local, pop.mutations,
                    offspring_data.second.mutation_keys, tables);

                // Give the caller a chance to generate
                // any metadata for the offspring that
                // may depend on the parents
                offspring_metadata[next_offspring].label = next_offspring;
                update_offspring(offspring_metadata[next_offspring], p1, p2,
                                 pop.diploid_metadata);
                // Update nodes of for offspring
                offspring_metadata[next_offspring].nodes[0]
                    = next_index_local - 1;
                offspring_metadata[next_offspring].nodes[1] = next_index_local;
            }
        assert(next_index_local
               == next_index + 2 * static_cast<std::int32_t>(N_next) - 1);
        // Copy the temp edges to the edge table
        auto x = temp_edges2.size();
        for (std::size_t i = 0; i < pop.diploids.size(); ++i)
            {
                temp_edges2.insert(end(temp_edges2), begin(temp_edges[i][0]),
                                   end(temp_edges[i][0]));
                temp_edges2.insert(end(temp_edges2), begin(temp_edges[i][1]),
                                   end(temp_edges[i][1]));
                temp_edges[i][0].clear();
                temp_edges[i][1].clear();
            }
        offsets.emplace_back(x, temp_edges2.size());
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
    }
} // namespace fwdpy11
#endif
