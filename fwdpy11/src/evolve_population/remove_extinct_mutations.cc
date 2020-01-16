#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <fwdpy11/types/Population.hpp>

namespace
{
    void
    reindex_container(const std::vector<fwdpp::uint_t>& indexes,
                      std::vector<fwdpp::uint_t>& keys)
    {
        for (auto& k : keys)
            {
                if (indexes[k] == std::numeric_limits<fwdpp::uint_t>::max())
                    {
                        throw std::runtime_error("bad mutation key remapping in haploid genome");
                    }
                k = indexes[k];
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
    std::vector<fwdpp::uint_t> new_mutation_indexes(
        pop.mutations.size(), std::numeric_limits<fwdpp::uint_t>::max());
    decltype(pop.mcounts) new_mcounts;
    decltype(pop.mcounts) new_preserved_mcounts;

    //Remove extinct mutations
    fwdpp::uint_t next_mutation_index = 0;
    for (fwdpp::uint_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (summed_counts[i] != 0)
                {
                    new_mutation_indexes[i] = next_mutation_index++;
                    new_mcounts.push_back(pop.mcounts[i]);
                    new_preserved_mcounts.push_back(
                        pop.mcounts_from_preserved_nodes[i]);
                }
            else
                {
                    pop.mutations[i].pos = std::numeric_limits<double>::max();
                }
        }
    pop.mutations.erase(
        std::remove_if(begin(pop.mutations), end(pop.mutations),
                       [](fwdpy11::Mutation& m) {
                           return m.pos == std::numeric_limits<double>::max();
                       }),
        end(pop.mutations));
    pop.mcounts.swap(new_mcounts);
    pop.mcounts_from_preserved_nodes.swap(new_preserved_mcounts);

    // Reindex the containers
    for (std::size_t i = 0; i < pop.tables.mutation_table.size(); ++i)
        {
            auto k = new_mutation_indexes[pop.tables.mutation_table[i].key];
            if (k == std::numeric_limits<fwdpp::uint_t>::max())
                {
                    throw std::runtime_error("bad mutation key remapping in mutation table");
                }
            pop.tables.mutation_table[i].key = k;
        }
    for (auto& g : pop.haploid_genomes)
        {
            if (g.n)
                {
                    reindex_container(new_mutation_indexes, g.mutations);
                    reindex_container(new_mutation_indexes, g.smutations);
                }
        }

    // Easiest way to update the lookup table:
    pop.mut_lookup.clear();
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            pop.mut_lookup.insert(std::make_pair(pop.mutations[i].pos, i));
        }
}

