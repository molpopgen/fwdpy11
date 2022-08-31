#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace
{
    template <typename GenomeVector>
    std::vector<std::size_t>
    compact_genomes(GenomeVector& input_genomes)
    {
        std::vector<std::size_t> genome_index(input_genomes.size(),
                                              std::numeric_limits<std::size_t>::max());
        std::size_t next_index = 0;
        for (std::size_t i = 0; i < input_genomes.size(); ++i)
            {
                if (input_genomes[i].n > 0)
                    {
                        genome_index[i] = next_index++;
                    }
            }
        GenomeVector genomes;
        genomes.reserve(next_index);
        for (std::size_t i = 0; i < input_genomes.size(); ++i)
            {
                if (genome_index[i] != std::numeric_limits<std::size_t>::max())
                    {
                        genomes.emplace_back(std::move(input_genomes[i]));
                    }
            }

        input_genomes.swap(genomes);
        return genome_index;
    }

    template <typename gcont_t>
    void
    validate(const std::size_t i, const gcont_t& haploid_genomes)
    {
        if (i == std::numeric_limits<std::size_t>::max() || i >= haploid_genomes.size())
            {
                throw std::runtime_error("error remapping genome indexes");
            }
        if (haploid_genomes[i].n == 0)
            {
                throw std::runtime_error("remapped genome is extinct");
            }
    }
} // namespace

void
remove_extinct_genomes(fwdpy11::DiploidPopulation& pop)
{
    auto genome_index = compact_genomes(pop.haploid_genomes);
    for (auto& dip : pop.diploids)
        {
            validate(genome_index[dip.first], pop.haploid_genomes);
            validate(genome_index[dip.second], pop.haploid_genomes);
            dip.first = genome_index[dip.first];
            dip.second = genome_index[dip.second];
        }
}

