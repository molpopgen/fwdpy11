#include <algorithm>
#include <numeric>
#include <fwdpy11/types/DiploidPopulation.hpp>

void
update(const fwdpp::haploid_genome &genome, std::vector<fwdpp::uint_t> &counts)
{
    for (auto k : genome.smutations)
        {
            counts[k]++;
        }
}

void
track_ancestral_counts(const std::vector<std::uint32_t> &individuals,
                       std::uint32_t *last_preserved_generation,
                       std::vector<std::uint32_t> &last_preserved_generation_counts,
                       fwdpy11::DiploidPopulation &pop)
{
    pop.mcounts_from_preserved_nodes.resize(pop.mutations.size(), 0);
    std::fill(begin(last_preserved_generation_counts),
              end(last_preserved_generation_counts), 0);
    last_preserved_generation_counts.resize(pop.mutations.size(), 0);
    for (auto i : individuals)
        {
            update(pop.haploid_genomes[pop.diploids[i].first],
                   last_preserved_generation_counts);
            update(pop.haploid_genomes[pop.diploids[i].second],
                   last_preserved_generation_counts);
        }
    *last_preserved_generation = pop.generation;
    std::transform(begin(last_preserved_generation_counts),
                   end(last_preserved_generation_counts),
                   begin(pop.mcounts_from_preserved_nodes),
                   begin(pop.mcounts_from_preserved_nodes), std::plus<std::uint32_t>());
}

