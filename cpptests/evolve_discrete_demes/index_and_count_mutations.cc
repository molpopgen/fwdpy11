#include <algorithm>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

void
index_and_count_mutations(bool suppress_edge_table_indexing,
                          bool simulating_neutral_variants,
                          bool reset_treeseqs_to_alive_nodes_after_simplification,
                          bool last_generation_was_recorded,
                          fwdpy11::DiploidPopulation& pop)
{
    pop.mcounts_from_preserved_nodes.resize(pop.mutations.size(), 0);
    if (!last_generation_was_recorded && !suppress_edge_table_indexing
        && !simulating_neutral_variants
        && !reset_treeseqs_to_alive_nodes_after_simplification)
        {
            return;
        }
    pop.tables->build_indexes();
    if (last_generation_was_recorded || pop.ancient_sample_metadata.empty()
        || simulating_neutral_variants)
        {
            pop.fill_alive_nodes();
            pop.fill_preserved_nodes();
            fwdpp::ts::count_mutations(*pop.tables, pop.mutations, pop.alive_nodes,
                                       pop.preserved_sample_nodes, pop.mcounts,
                                       pop.mcounts_from_preserved_nodes);
        }
    else
        {
            fwdpp::fwdpp_internal::process_haploid_genomes(pop.haploid_genomes,
                                                           pop.mutations, pop.mcounts);
        }
    if (reset_treeseqs_to_alive_nodes_after_simplification)
        {
            auto itr = std::remove_if(
                begin(pop.tables->mutations), end(pop.tables->mutations),
                [&pop](const auto& mr) {
                    return pop.mcounts[mr.key] + pop.mcounts_from_preserved_nodes[mr.key]
                           == 0;
                });
            auto d = std::distance(itr, end(pop.tables->mutations));
            pop.tables->mutations.erase(itr, end(pop.tables->mutations));
            if (d)
                {
                    fwdpp::ts::rebuild_site_table(*pop.tables);
                }
        }
}

