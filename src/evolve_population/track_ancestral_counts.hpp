#ifndef FWDPY11_TSEVOLVE_TRACK_ANCESTRAL_COUNTS_HPP
#define FWDPY11_TSEVOLVE_TRACK_ANCESTRAL_COUNTS_HPP

#include <cstdint>
#include <vector>

namespace fwdpy11
{
    class DiploidPopulation;
}

void track_ancestral_counts(const std::vector<std::uint32_t> &individuals,
                            std::uint32_t *last_preserved_generation,
                            std::vector<std::uint32_t> &last_preserved_generation_counts,
                            fwdpy11::DiploidPopulation &pop);

#endif
