#include <vector>
#include <algorithm>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpy11/types/Diploid.hpp>

void
cleanup_metadata(const fwdpp::ts::table_collection& tables,
                 const fwdpp::uint_t generation,
                 std::vector<fwdpy11::DiploidMetadata>& metadata)
{
    metadata.erase(
        std::remove_if(begin(metadata), end(metadata),
                       [generation, &tables](fwdpy11::DiploidMetadata& md) {
                           return tables.node_table[md.nodes[0]].time
                                      == generation
                                  || tables.node_table[md.nodes[1]].time
                                         == generation;
                       }),
        end(metadata));
}

