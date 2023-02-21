#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/ts/definitions.hpp>

void
remap_metadata(std::vector<fwdpy11::DiploidMetadata> &metadata,
               const std::vector<fwdpp::ts::table_index_t> &idmap)
{
    for (auto &m : metadata)
        {
            m.nodes[0] = idmap[m.nodes[0]];
            m.nodes[1] = idmap[m.nodes[1]];
            if (m.nodes[0] == fwdpp::ts::NULL_INDEX
                || m.nodes[1] == fwdpp::ts::NULL_INDEX)
                {
                    throw std::runtime_error(
                        "error remapping node field of individual metadata");
                }
        }
}

void
coordinate_count_vector_sizes(
    std::size_t nmutations, std::vector<std::uint32_t> &mcounts,
    std::vector<std::uint32_t> &mcounts_from_preserved_nodes)
{
    mcounts.resize(nmutations, 0);
    mcounts_from_preserved_nodes.resize(nmutations, 0);
}
