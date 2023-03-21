#pragma once

#include <fwdpy11/types/DiploidPopulation.hpp>

/* This is the back end for importing mutations from tskit.
 * They have been decoded from metadata already.
 *
 * The procedure is:
 * * Add them to our table structures.
 * * Sort the tables.
 * * Traverse trees so that we can add mutations to genomes.
 *
 * Because the data come from a tskit TreeSequence, we
 * can make some simplifying assumptions:
 *
 * 1. Edges, etc., are already sorted.
 * 2. So we just need to sort mutations.
 */
void
set_mutations(const std::vector<fwdpy11::Mutation> &mutations,
              const std::vector<std::int32_t> &mutation_nodes,
              const std::vector<fwdpy11::mutation_origin_time> &origin_times,
              fwdpy11::DiploidPopulation &pop);
