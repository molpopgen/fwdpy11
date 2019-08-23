#ifndef FWDPY11_TSEVOLVE_TRACK_ANCESTRAL_COUNTS_HPP
#define FWDPY11_TSEVOLVE_TRACK_ANCESTRAL_COUNTS_HPP

#include <fwdpy11/types/DiploidPopulation.hpp>

void track_ancestral_counts(fwdpy11::DiploidPopulation &pop,
                            const std::vector<fwdpp::uint_t> &individuals);

#endif
