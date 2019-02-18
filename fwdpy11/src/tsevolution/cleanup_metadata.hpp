#ifndef FWDPY11_TSEVOLUTION_CLEANUP_METADATA_HPP
#define FWDPY11_TSEVOLUTION_CLEANUP_METADATA_HPP

#include <vector>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpy11/types/Diploid.hpp>

void
cleanup_metadata(const fwdpp::ts::table_collection& tables,
                 const fwdpp::uint_t generation,
                 std::vector<fwdpy11::DiploidMetadata>& metadata);

#endif
