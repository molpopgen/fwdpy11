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
track_ancestral_counts(fwdpy11::DiploidPopulation &pop,
                       const std::vector<fwdpp::uint_t> &individuals)
{
    pop.mcounts_from_preserved_nodes.resize(pop.mutations.size(), 0);
    for (auto i : individuals)
        {
            update(pop.haploid_genomes[pop.diploids[i].first],
                   pop.mcounts_from_preserved_nodes);
            update(pop.haploid_genomes[pop.diploids[i].second],
                   pop.mcounts_from_preserved_nodes);
        }
}

