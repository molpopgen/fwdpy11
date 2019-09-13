#include <fwdpy11/evolvets/simplify_tables.hpp>
#include "reset_treeseq_to_alive_nodes.hpp"
#include "util.hpp"

void
reset_treeseqs_to_alive_nodes(
    const fwdpy11::DiploidPopulation_temporal_sampler& recorder,
    const bool preserve_selected_fixations,
    const bool suppress_edge_table_indexing,
    fwdpy11::DiploidPopulation& pop,
    fwdpp::ts::table_simplifier&simplifier,
    fwdpp::flagged_mutation_queue &mutation_recycling_bin)
{
    recorder(pop);
    if (!pop.tables.preserved_nodes.empty())
        {
            pop.tables.preserved_nodes.clear();
            pop.ancient_sample_metadata.clear();
            pop.ancient_sample_records.clear();
            pop.ancient_sample_genetic_value_matrix.clear();
            auto rv = fwdpy11::simplify_tables(
                pop, pop.mcounts_from_preserved_nodes, pop.tables, simplifier,
                preserve_selected_fixations, false,
                suppress_edge_table_indexing);
            if (suppress_edge_table_indexing == false)
                {
                    mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                        pop.mcounts, pop.mcounts_from_preserved_nodes);
                }
            else
                {
                    mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                        rv.second, pop.mutations.size());
                }
            remap_metadata(pop.ancient_sample_metadata, rv.first);
            remap_metadata(pop.diploid_metadata, rv.first);
        }
}
