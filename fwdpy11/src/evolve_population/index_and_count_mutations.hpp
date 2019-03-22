#ifndef FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP
#define FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP

#include <vector>
#include <cstdint>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpp/ts/table_collection.hpp>

void
index_and_count_mutations(
    const bool suppress_edge_table_indexing, const std::size_t twoN,
    const std::vector<fwdpy11::Mutation>& mutations,
    fwdpp::ts::table_collection& tables,
    std::vector<fwdpp::uint_t>& mcounts,
    std::vector<fwdpp::uint_t>& mcounts_from_preserved_nodes);

#endif
