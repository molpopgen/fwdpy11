#ifndef FWDPY11_TSEVOLVE_UTIL_HPP
#define FWDPY11_TSEVOLVE_UTIL_HPP

#include <cstdint>
#include <vector>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/ts/definitions.hpp>

void
remap_metadata(std::vector<fwdpy11::DiploidMetadata> &metadata,
               const std::vector<fwdpp::ts::TS_NODE_INT> &idmap);
void
coordinate_count_vector_sizes(
    std::size_t nmutations, std::vector<std::uint32_t> &mcounts,
    std::vector<std::uint32_t> &mcounts_from_preserved_nodes);

#endif
