#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>
#include <fwdpy11/types/Population.hpp>

namespace
{
    std::size_t
    reindex(const std::vector<std::size_t>& indexes,
            const std::size_t input_index)
    {
        auto f = std::find(begin(indexes), end(indexes), input_index);
        auto d = std::distance(begin(indexes), f);
        return d;
    }

    void
    reindex_container(const std::vector<std::size_t>& indexes,
                      std::vector<fwdpp::uint_t>& keys)
    {
        for (auto& k : keys)
            {
                k = reindex(indexes, k);
            }
    }
} // namespace

void
remove_extinct_mutations(fwdpy11::Population& pop)
{
    std::vector<fwdpp::uint_t> summed_counts(pop.mcounts);
    std::transform(begin(pop.mcounts_from_preserved_nodes),
                   end(pop.mcounts_from_preserved_nodes), begin(summed_counts),
                   begin(summed_counts), std::plus<fwdpp::uint_t>());
    std::vector<std::size_t> indexes(pop.mutations.size());
    std::iota(begin(indexes), end(indexes), 0);

    //Remove extinct mutations
    for (auto i : indexes)
        {
            if (summed_counts[i] == 0)
                {
                    pop.mutations[i].pos = std::numeric_limits<double>::max();
                    pop.mcounts[i] = std::numeric_limits<fwdpp::uint_t>::max();
                    pop.mcounts_from_preserved_nodes[i]
                        = std::numeric_limits<fwdpp::uint_t>::max();
                }
        }
    pop.mutations.erase(
        std::remove_if(begin(pop.mutations), end(pop.mutations),
                       [](fwdpy11::Mutation& m) {
                           return m.pos == std::numeric_limits<double>::max();
                       }),
        end(pop.mutations));
    pop.mcounts.erase(std::remove(begin(pop.mcounts), end(pop.mcounts),
                                  std::numeric_limits<fwdpp::uint_t>::max()),
                      end(pop.mcounts));
    pop.mcounts_from_preserved_nodes.erase(
        std::remove(begin(pop.mcounts_from_preserved_nodes),
                    end(pop.mcounts_from_preserved_nodes),
                    std::numeric_limits<fwdpp::uint_t>::max()),
        end(pop.mcounts_from_preserved_nodes));

    // Erase indexes associated w/extinct mutations
    indexes.erase(std::remove_if(begin(indexes), end(indexes),
                                 [&summed_counts](const std::size_t i) {
                                     return summed_counts[i] == 0;
                                 }),
                  end(indexes));

    // Reindex the containers
    for (std::size_t i = 0; i < pop.tables.mutation_table.size(); ++i)
        {
            pop.tables.mutation_table[i].key
                = reindex(indexes, pop.tables.mutation_table[i].key);
        }
    for (auto& g : pop.gametes)
        {
            if (g.n)
                {
                    reindex_container(indexes, g.mutations);
                    reindex_container(indexes, g.smutations);
                }
        }

    // Easiest way to update the lookup table:
    pop.mut_lookup.clear();
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            pop.mut_lookup.insert(std::make_pair(pop.mutations[i].pos, i));
        }
}

