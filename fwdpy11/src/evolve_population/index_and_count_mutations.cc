#include <numeric>
#include <fwdpp/ts/count_mutations.hpp>
#include "index_and_count_mutations.hpp"

void
index_and_count_mutations(
    const bool suppress_edge_table_indexing, const std::size_t twoN,
    const std::vector<fwdpy11::Mutation>& mutations,
    fwdpp::ts::table_collection& tables,
    std::vector<fwdpp::uint_t>& mcounts,
    std::vector<fwdpp::uint_t>& mcounts_from_preserved_nodes)
{
    if (!suppress_edge_table_indexing)
        {
            return;
        }
    tables.build_indexes();
    std::vector<std::int32_t> samples(twoN);
    std::iota(samples.begin(), samples.end(), 0);
    fwdpp::ts::count_mutations(tables, mutations, samples, mcounts,
                               mcounts_from_preserved_nodes);
}

