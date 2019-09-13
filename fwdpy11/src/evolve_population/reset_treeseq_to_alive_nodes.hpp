#ifndef FWDPY11_TSEVOLUTION_RESET_TREESEQS_TO_ALIVE_NODES_HPP
#define FWDPY11_TSEVOLUTION_RESET_TREESEQS_TO_ALIVE_NODES_HPP

#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/ts/table_simplifier.hpp>

void reset_treeseqs_to_alive_nodes(
    const fwdpy11::DiploidPopulation_temporal_sampler& recorder,
    const bool preserve_selected_fixations,
    const bool suppress_edge_table_indexing,
    fwdpy11::DiploidPopulation& pop,
    fwdpp::ts::table_simplifier&simplifier,
    fwdpp::flagged_mutation_queue &mutation_recycling_bin 
    );

#endif

