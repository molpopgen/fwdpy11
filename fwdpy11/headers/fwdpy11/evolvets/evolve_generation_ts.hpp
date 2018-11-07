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
// This file is copied from the fwdpp version 0.7.0 examples folder.
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
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>

namespace fwdpy11
{

    // Wow, that's a lot of stuff needed:
    template <typename poptype, typename rng_t, typename breakpoint_function,
              typename mutation_model, typename mrecbin, typename grecbin>
    std::int32_t
    generate_offspring(
        const rng_t& rng, const breakpoint_function& recmodel,
        const mutation_model& mmodel, const double mu,
        const std::size_t parent, const fwdpp::uint_t parent_g1,
        const fwdpp::uint_t parent_g2,
        const std::tuple<std::int32_t, std::int32_t>& parent_nodes,
        const std::int32_t generation, const std::int32_t next_index,
        poptype& pop, std::size_t& offspring_gamete,
        fwdpp::ts::table_collection& tables, mrecbin& mutation_recycling_bin,
        grecbin& gamete_recycling_bin)
    {
        auto breakpoints = recmodel();
        auto new_mutations = fwdpp::generate_new_mutations(
            mutation_recycling_bin, rng.get(), mu, pop.diploids[parent],
            pop.gametes, pop.mutations, parent_g1, mmodel);
#ifndef NDEBUG
        for (auto& m : new_mutations)
            {
                auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
                assert(std::distance(itr.first, itr.second) == 1);
            }
#endif
        // We will only add selected mutations into offspring gametes.
        auto end_of_neutral = std::stable_partition(
            new_mutations.begin(), new_mutations.end(),
            [&pop](const fwdpp::uint_t key) {
                return pop.mutations[key].neutral == true;
            });
        offspring_gamete = fwdpp::mutate_recombine(
            decltype(new_mutations)(end_of_neutral, new_mutations.end()),
            breakpoints, parent_g1, parent_g2, pop.gametes, pop.mutations,
            gamete_recycling_bin, pop.neutral, pop.selected);
        //TODO: generalize this for offspring deme != 0
        tables.add_offspring_data(next_index, breakpoints, new_mutations,
                                  parent_nodes, 0, generation);
        return next_index + 1;
    }

    //TODO: need to track contribution to ancestral mutation counts!
    template <typename rng_t, typename poptype, typename pick_parent1_fxn,
              typename pick_parent2_fxn, typename offspring_metadata_fxn,
              typename breakpoint_function, typename mutation_model>
    void
    evolve_generation_ts(const rng_t& rng, poptype& pop,
                         const fwdpp::uint_t N_next, const double mu,
                         const pick_parent1_fxn& pick1,
                         const pick_parent2_fxn& pick2,
                         const offspring_metadata_fxn& update_offspring,
                         const mutation_model& mmodel,
                         std::queue<std::size_t>& mutation_recycling_bin,
                         const breakpoint_function& recmodel,
                         const fwdpp::uint_t generation,
                         fwdpp::ts::table_collection& tables,
                         std::int32_t first_parental_index,
                         std::int32_t next_index)
    {

        auto gamete_recycling_bin
            = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
        for (auto& dip : pop.diploids)
            {
                pop.gametes[dip.first].n = pop.gametes[dip.second].n = 0;
            }

        decltype(pop.diploids) offspring(N_next);
        decltype(pop.diploid_metadata) offspring_metadata(N_next);

        // Generate the offspring
        auto next_index_local = next_index;
        for (std::size_t next_offspring = 0; next_offspring < offspring.size();
             ++next_offspring)
            {
                auto p1 = pick1();
                auto p2 = pick2(p1);
                auto p1g1 = pop.diploids[p1].first;
                auto p1g2 = pop.diploids[p1].second;
                auto p2g1 = pop.diploids[p2].first;
                auto p2g2 = pop.diploids[p2].second;

                // Mendel
                int swap1 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
                int swap2 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
                if (swap1)
                    std::swap(p1g1, p1g2);
                if (swap2)
                    std::swap(p2g1, p2g2);

                auto p1id = fwdpp::ts::get_parent_ids(first_parental_index, p1,
                                                      swap1);
                auto p2id = fwdpp::ts::get_parent_ids(first_parental_index, p2,
                                                      swap2);

                auto& dip = offspring[next_offspring];
                next_index_local = generate_offspring(
                    rng, recmodel, mmodel, mu, p1, p1g1, p1g2, p1id,
                    generation, next_index_local, pop, dip.first, tables,
                    mutation_recycling_bin, gamete_recycling_bin);
                next_index_local = generate_offspring(
                    rng, recmodel, mmodel, mu, p2, p2g1, p2g2, p2id,
                    generation, next_index_local, pop, dip.second, tables,
                    mutation_recycling_bin, gamete_recycling_bin);
                pop.gametes[dip.first].n++;
                pop.gametes[dip.second].n++;
                // Give the caller a chance to generate
                // any metadata for the offspring that
                // may depend on the parents
                offspring_metadata[next_offspring].label = next_offspring;
                update_offspring(offspring_metadata[next_offspring], p1, p2,
                                 pop.diploid_metadata);
            }
        assert(next_index_local
               == next_index + 2 * static_cast<std::int32_t>(N_next));
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
    }
} // namespace fwdpy11
#endif
