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
#include <fwdpp/ts/recycling.hpp>
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
                    const fwdpp::ts::TS_NODE_INT first_sample_node,
                    const std::size_t num_samples,
                    const bool preserve_selected_fixations,
                    const bool simulating_neutral_variants,
                    const bool suppress_edge_table_indexing)
    {
        tables.sort_tables(pop.mutations);
        std::vector<std::int32_t> samples(num_samples);
        std::iota(samples.begin(), samples.end(), first_sample_node);
        auto rv = simplifier.simplify(tables, samples, pop.mutations);

        for (auto &s : samples)
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
                return rv;
            }
        tables.build_indexes();
        fwdpp::ts::count_mutations(tables, pop.mutations, samples, pop.mcounts,
                                   mcounts_from_preserved_nodes);
        // TODO: update this to allow neutral mutations to be simulated
        // TODO: better fixation handling via accounting for number of ancient samples
        if (!preserve_selected_fixations && !simulating_neutral_variants)
            {
                tables.mutation_table.erase(
                    std::remove_if(
                        tables.mutation_table.begin(),
                        tables.mutation_table.end(),
                        [&pop, &mcounts_from_preserved_nodes](
                            const fwdpp::ts::mutation_record &mr) {
                            return pop.mcounts[mr.key]
                                       == 2 * pop.diploids.size()
                                   && mcounts_from_preserved_nodes[mr.key]
                                          == 0;
                        }),
                    tables.mutation_table.end());
                fwdpp::ts::remove_fixations_from_gametes(
                    pop.gametes, pop.mutations, pop.mcounts,
                    mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                    preserve_selected_fixations);
            }

        // TODO: the blocks below should be abstracted out into a closure
        // that removes these if statements
        if (!preserve_selected_fixations && !simulating_neutral_variants)
            {
                fwdpp::ts::flag_mutations_for_recycling(
                    pop, mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                    pop.generation, std::false_type(), std::false_type());
            }
        else if (preserve_selected_fixations && !simulating_neutral_variants)
            {
                fwdpp::ts::flag_mutations_for_recycling(
                    pop, mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                    pop.generation, std::true_type(), std::false_type());
            }
        //confirm_mutation_counts(pop, tables);
        return rv;
    }
} // namespace fwdpy11
#endif
