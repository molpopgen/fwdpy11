#ifndef FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP
#define FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP

#include <vector>
#include <cstdint>
#include <fwdpy11/types/DiploidPopulation.hpp>

void index_and_count_mutations(const bool suppress_edge_table_indexing,
                               fwdpy11::DiploidPopulation& pop);

#endif
