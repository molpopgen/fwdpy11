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

    template <typename rng_t, typename poptype, typename pick_parent1_fxn,
              typename pick_parent2_fxn, typename offspring_metadata_fxn,
              typename genetic_param_holder>
    void
    evolve_generation_ts(
        const rng_t& rng, poptype& pop, genetic_param_holder& genetics,
        const fwdpp::uint_t N_next, const pick_parent1_fxn& pick1,
        const pick_parent2_fxn& pick2,
        const offspring_metadata_fxn& update_offspring,
        const fwdpp::uint_t generation, fwdpp::ts::table_collection& tables,
        std::int32_t first_parental_index, std::int32_t next_index)
    {
        fwdpp::debug::all_gametes_extant(pop);

        genetics.gamete_recycling_bin = fwdpp::make_gamete_queue(pop.gametes);

        fwdpp::zero_out_gametes(pop);

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
                tables.add_offspring_data(
                    next_index_local++, offspring_data.first.breakpoints,
                    offspring_data.first.mutation_keys, p1id, 0, generation);
                tables.add_offspring_data(
                    next_index_local++, offspring_data.second.breakpoints,
                    offspring_data.second.mutation_keys, p2id, 0, generation);

                // Give the caller a chance to generate
                // any metadata for the offspring that
                // may depend on the parents
                offspring_metadata[next_offspring].label = next_offspring;
                update_offspring(offspring_metadata[next_offspring], p1, p2,
                                 pop.diploid_metadata);
                // Update nodes of for offspring
                offspring_metadata[next_offspring].nodes[0]
                    = next_index_local - 2;
                offspring_metadata[next_offspring].nodes[1]
                    = next_index_local - 1;
            }
        assert(next_index_local
               == next_index + 2 * static_cast<std::int32_t>(N_next));
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
    }
} // namespace fwdpy11
#endif
