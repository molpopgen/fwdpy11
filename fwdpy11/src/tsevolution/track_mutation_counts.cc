#include <tuple>
#include <algorithm>
#include <fwdpy11/types/Population.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

void
track_mutation_counts(fwdpy11::Population &pop, const bool simplified,
                      const bool suppress_edge_table_indexing)
{
    if (!simplified || (simplified && suppress_edge_table_indexing))
        {
            fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                                   pop.mcounts);
        }
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
