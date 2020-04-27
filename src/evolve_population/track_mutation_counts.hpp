#ifndef FWDPY11_TSEVOLUTION_TRACK_MUTATION_COUNTS_HPP
#define FWDPY11_TSEVOLUTION_TRACK_MUTATION_COUNTS_HPP

#include <fwdpy11/types/Population.hpp>

void track_mutation_counts(fwdpy11::Population &pop,
                           const bool simplified,
                           const bool suppress_edge_table_indexing);

#endif
