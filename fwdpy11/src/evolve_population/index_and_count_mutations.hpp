#ifndef FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP
#define FWDPP_TSEVOLUTION_INDEX_AND_COUNT_HPP

#include <fwdpy11/types/DiploidPopulation.hpp>

void index_and_count_mutations(bool suppress_edge_table_indexing,
                               bool simulating_neutral_variants,
                               bool reset_treeseqs_to_alive_nodes_after_simplification,
                               fwdpy11::DiploidPopulation& pop);

#endif
