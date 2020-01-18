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
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpy11/discrete_demography/simulation.hpp>

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

    inline std::pair<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>
    parent_nodes_from_metadata(
        const std::size_t i,
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
        // NOTE: could the manager? be const?
        fwdpy11::discrete_demography::discrete_demography_manager&
            ddemog_manager,
        const fwdpp::uint_t generation, fwdpp::ts::table_collection& tables,
        std::int32_t next_index)
    {
        fwdpp::debug::all_haploid_genomes_extant(pop);

        genetics.haploid_genome_recycling_bin
            = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);

        fwdpp::zero_out_haploid_genomes(pop);

        decltype(pop.diploids) offspring;
        decltype(pop.diploid_metadata) offspring_metadata;
        offspring.reserve(pop.N);
        offspring_metadata.reserve(pop.N);

        // Generate the offspring
        auto next_index_local = next_index;
        using maxdeme_type = typename std::remove_const<decltype(
            ddemog_manager.maxdemes)>::type;
        for (maxdeme_type deme = 0; deme < ddemog_manager.maxdemes; ++deme)
            {
                auto next_N_deme
                    = ddemog_manager.sizes_rates.next_deme_sizes.get()[deme];
                for (decltype(next_N_deme) ind = 0; ind < next_N_deme; ++ind)
                    {
                        // Get the parents
                        auto pdata
                            = fwdpy11::discrete_demography::pick_parents(
                                rng, deme, ddemog_manager.miglookup,
                                ddemog_manager.sizes_rates.current_deme_sizes,
                                ddemog_manager.sizes_rates.selfing_rates,
                                ddemog_manager.fitnesses);
                        // Add a new diploid
                        offspring.emplace_back(fwdpy11::DiploidGenotype{
                            std::numeric_limits<std::size_t>::max(),
                            std::numeric_limits<std::size_t>::max() });
                        auto& dip = offspring.back();
                        //auto p1 = pick1();
                        //auto p2 = pick2(p1);
                        auto offspring_data = generate_offspring(
                            rng, std::make_pair(pdata.parent1, pdata.parent2),
                            pop, dip, genetics);
                        auto p1id = parent_nodes_from_metadata(
                            pdata.parent1, pop.diploid_metadata,
                            offspring_data.first.swapped);
                        auto p2id = parent_nodes_from_metadata(
                            pdata.parent2, pop.diploid_metadata,
                            offspring_data.second.swapped);
                        fwdpp::ts::TS_NODE_INT offspring_node_1
                            = tables.register_diploid_offspring(
                                offspring_data.first.breakpoints, p1id, deme,
                                generation);
                        fwdpp::ts::record_mutations_infinite_sites(
                            offspring_node_1, pop.mutations,
                            offspring_data.first.mutation_keys, tables);
                        fwdpp::ts::TS_NODE_INT offspring_node_2
                            = tables.register_diploid_offspring(
                                offspring_data.second.breakpoints, p2id, deme,
                                generation);
                        fwdpp::ts::record_mutations_infinite_sites(
                            offspring_node_2, pop.mutations,
                            offspring_data.second.mutation_keys, tables);

                        // Add metadata for the offspring
                        offspring_metadata.emplace_back(
                            fwdpy11::DiploidMetadata{
                                0.0,
                                0.0,
                                1.,
                                { 0, 0, 0 },
                                offspring_metadata.size(),
                                { pdata.parent1, pdata.parent2 },
                                deme,
                                0,
                                { offspring_node_1, offspring_node_2 } });

                        next_index_local = offspring_node_2;
                    }
            }
        assert(next_index_local == pop.tables.num_nodes() - 1);
        if (next_index_local
            != static_cast<decltype(next_index_local)>(pop.tables.num_nodes()
                                                       - 1))
            {
                throw std::runtime_error(
                    "error in book-keeping offspring nodes");
            }
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
        pop.N = static_cast<std::uint32_t>(pop.diploid_metadata.size());
    }
} // namespace fwdpy11
#endif
