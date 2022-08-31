#include <tuple>
#include <algorithm>
#include <fwdpy11/types/Population.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include "util.hpp"

void
track_mutation_counts(fwdpy11::Population &pop, const bool simplified,
                      const bool suppress_edge_table_indexing)
{
    if (pop.mcounts.size() != pop.mcounts_from_preserved_nodes.size())
        {
            throw std::runtime_error(
                "track_mutation_counts: count vector size mismatch");
        }
    if (!simplified || suppress_edge_table_indexing)
        {
            fwdpp::fwdpp_internal::process_haploid_genomes(
                pop.haploid_genomes, pop.mutations, pop.mcounts);
        }
    coordinate_count_vector_sizes(pop.mutations.size(), pop.mcounts,
                                  pop.mcounts_from_preserved_nodes);
    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
        {
            if (pop.mcounts[i] == 2 * pop.N)
                {
                    auto loc = std::lower_bound(
                        pop.fixations.begin(), pop.fixations.end(),
                        std::make_tuple(pop.mutations[i].g,
                                        pop.mutations[i].pos),
                        [](const fwdpy11::Mutation &mut,
                           const std::tuple<double, std::uint32_t>
                               &value) noexcept {
                            return std::tie(mut.g, mut.pos) < value;
                        });
                    auto d = std::distance(pop.fixations.begin(), loc);
                    if (loc == end(pop.fixations)
                        || (loc->pos != pop.mutations[i].pos
                            && loc->g != pop.mutations[i].g))
                        {
                            pop.fixations.insert(loc, pop.mutations[i]);
                            pop.fixation_times.insert(
                                begin(pop.fixation_times) + d, pop.generation);
                        }
                }
        }
}
