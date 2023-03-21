//
// Copyright (C) 2017-2023 Kevin Thornton <krthornt@uci.edu>
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

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/recording/edge_buffer.hpp>
#include <fwdpp/ts/simplify_tables.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

// TODO allow for fixation recording
// and simulation of neutral variants
template <typename poptype, typename SimplificationState>
void
simplify_tables(poptype &pop, std::vector<fwdpp::uint_t> &mcounts_from_preserved_nodes,
                std::vector<fwdpp::ts::table_index_t> &alive_at_last_simplification,
                fwdpp::ts::std_table_collection &tables,
                SimplificationState &simplifier_state,
                fwdpp::ts::simplify_tables_output &simplification_output,
                fwdpp::ts::edge_buffer &new_edge_buffer,
                const bool preserve_selected_fixations,
                const bool suppress_edge_table_indexing)
{
    // As of 0.8.0, we do not need to sort edges!
    fwdpp::ts::sort_mutation_table(tables);
    pop.fill_alive_nodes();
    pop.fill_preserved_nodes();
    auto samples(pop.alive_nodes);
    samples.insert(end(samples), begin(pop.preserved_sample_nodes),
                   end(pop.preserved_sample_nodes));
    fwdpp::ts::simplify_tables(samples, alive_at_last_simplification,
                               fwdpp::ts::simplification_flags{}, simplifier_state,
                               *pop.tables, new_edge_buffer, simplification_output);
    for (auto &s : pop.alive_nodes)
        {
            s = simplification_output.idmap[s];
        }
    for (auto &s : pop.preserved_sample_nodes)
        {
            s = simplification_output.idmap[s];
        }
    // Remove mutations that are simplified out
    // from the population hash table.
    std::vector<int> preserved(pop.mutations.size(), 0);
    for (auto p : simplification_output.preserved_mutations)
        {
            preserved[p] = 1;
        }
    for (std::size_t p = 0; p < preserved.size(); ++p)
        {
            if (!preserved[p])
                {
                    fwdpp::ts::detail::process_mutation_index(pop.mutations,
                                                              pop.mut_lookup, p);
                }
        }

    if (suppress_edge_table_indexing == true)
        {
            pop.mcounts.resize(pop.mutations.size(), 0);
            pop.mcounts_from_preserved_nodes.resize(pop.mutations.size(), 0);
            return;
        }
    pop.tables->build_indexes();
    if (pop.ancient_sample_metadata.empty())
        {
            fwdpp::ts::count_mutations(tables, pop.mutations, pop.alive_nodes,
                                       pop.preserved_sample_nodes, pop.mcounts,
                                       mcounts_from_preserved_nodes);
        }
    else
        {
            fwdpp::fwdpp_internal::process_haploid_genomes(pop.haploid_genomes,
                                                           pop.mutations, pop.mcounts);
            pop.mcounts_from_preserved_nodes.resize(pop.mcounts.size(), 0);
        }

    // TODO: revisit this comment.
    // in pr 779, we removed simulating_neutral_variants b/c it was unused.
    //
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
        pop.tables->mutations.begin(), pop.tables->mutations.end(),
        [&pop, &mcounts_from_preserved_nodes,
         preserve_selected_fixations](const fwdpp::ts::mutation_record &mr) {
            if (pop.mutations[mr.key].neutral == false && preserve_selected_fixations)
                {
                    return false;
                }
            return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                   && mcounts_from_preserved_nodes[mr.key] == 0;
        });
    auto d = std::distance(itr, end(pop.tables->mutations));
    pop.tables->mutations.erase(itr, end(pop.tables->mutations));
    if (d)
        {
            fwdpp::ts::rebuild_site_table(tables);
        }
}
