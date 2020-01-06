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
// Copied from fwdpp 0.7.0 examples.  Author is Kevin Thornton
#ifndef FWDPY11_SIMPLIFY_TABLES_HPP
#define FWDPY11_SIMPLIFY_TABLES_HPP

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
//#include "confirm_mutation_counts.hpp"

namespace fwdpy11
{

    // TODO allow for fixation recording
    // and simulation of neutral variants
    template <typename poptype>
    std::pair<std::vector<fwdpp::ts::TS_NODE_INT>, std::vector<std::size_t>>
    simplify_tables(poptype &pop,
                    std::vector<fwdpp::uint_t> &mcounts_from_preserved_nodes,
                    fwdpp::ts::table_collection &tables,
                    fwdpp::ts::table_simplifier &simplifier,
                    const bool preserve_selected_fixations,
                    const bool simulating_neutral_variants,
                    const bool suppress_edge_table_indexing)
    {
        tables.sort_tables_for_simplification();
        pop.fill_alive_nodes();
        auto rv = simplifier.simplify(tables, pop.alive_nodes);

        for (auto &s : pop.alive_nodes)
            {
                s = rv.first[s];
            }
#ifndef NDEBUG
        for (auto &s : tables.preserved_nodes)
            {
                if (s == -1)
                    {
                        throw std::runtime_error("ancient sample node is NULL "
                                                 "after simplification");
                    }
            }
#endif
        if (suppress_edge_table_indexing == true)
            {
                pop.mcounts.resize(pop.mutations.size(), 0);
                pop.mcounts_from_preserved_nodes.resize(pop.mutations.size(), 0);
                return rv;
            }
        tables.build_indexes();
        if (pop.tables.preserved_nodes.empty())
            {
                fwdpp::ts::count_mutations(tables, pop.mutations,
                                           pop.alive_nodes, pop.mcounts,
                                           mcounts_from_preserved_nodes);
            }
        else
            {
                fwdpp::fwdpp_internal::process_haploid_genomes(
                    pop.haploid_genomes, pop.mutations, pop.mcounts);
                pop.mcounts_from_preserved_nodes.resize(pop.mcounts.size(), 0);
            }

        // If we are here, then the tables are indexed and mutations are counted.
        // If there are ancient samples and we are simulating neutral variants,
        // then we do not know their frequencies.  Thus, neutral
        // fixations will not be removed until the end of the simulation.
        // To change this behavior, the logic above would need to allow
        // for mutation counting from tree seqs w/ancient samples or to
        // force a tree traversal here if there are ancient samples.
        // Such a traversal would allow for the purging of "global" fixations,
        // e.g., mutations on root nodes in trees w/only 1 root.
        // NOTE: if we use tree sequences to globablly purge fixations,
        // then we cannot test that genetic values of ancient sample metadata
        // equal genetic values calculated from tree sequences.
        // Behavior change in 0.5.3: check for fixations
        // no matter what.

        auto itr = std::remove_if(
            tables.mutation_table.begin(), tables.mutation_table.end(),
            [&pop, &mcounts_from_preserved_nodes, preserve_selected_fixations](
                const fwdpp::ts::mutation_record &mr) {
                if (pop.mutations[mr.key].neutral == false
                    && preserve_selected_fixations)
                    {
                        return false;
                    }
                return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                       && mcounts_from_preserved_nodes[mr.key] == 0;
            });
        auto d = std::distance(itr, end(tables.mutation_table));
        tables.mutation_table.erase(itr, end(tables.mutation_table));
        if (d)
            {
                tables.rebuild_site_table();
            }
        return rv;
    } // namespace fwdpy11
} // namespace fwdpy11
#endif
