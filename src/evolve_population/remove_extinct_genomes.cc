#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace
{
    template <typename gcont_t>
    std::size_t
    process_genome(const std::size_t current_index, gcont_t& input_genomes,
                   gcont_t& output_genomes,
                   std::vector<std::size_t>& genome_index,
                   std::vector<int>& processed)
    {
        std::size_t rv = genome_index[current_index];
        if (!processed[current_index])
            {
                output_genomes.emplace_back(
                    std::move(input_genomes[current_index]));
                processed[current_index] = 1;
                genome_index[current_index] = output_genomes.size() - 1;
                rv = genome_index[current_index];
            }
        return rv;
    }

    template<typename gcont_t> 
    void validate(const std::size_t i, const gcont_t & haploid_genomes)
    {
        if(i ==std::numeric_limits<std::size_t>::max() || i >= haploid_genomes.size())
        {
            throw std::runtime_error("error remapping genome indexes");
        }
        if(haploid_genomes[i].n == 0)
        {
            throw std::runtime_error("remapped genome is extinct");
        }
    }
} // namespace

void
remove_extinct_genomes(fwdpy11::DiploidPopulation& pop)
{
    decltype(pop.haploid_genomes) genomes;
    std::vector<std::size_t> genome_index(
        pop.haploid_genomes.size(), std::numeric_limits<std::size_t>::max());
    std::vector<int> processed(genome_index.size(), 0);
    std::vector<int> isizes(genome_index.size(),0);
    for (auto& dip : pop.diploids)
        {
            dip.first = process_genome(dip.first, pop.haploid_genomes, genomes,
                                       genome_index, processed);
            validate(dip.first, genomes);
            dip.second = process_genome(dip.second, pop.haploid_genomes,
                                        genomes, genome_index, processed);
            validate(dip.second, genomes);
        }
    pop.haploid_genomes.swap(genomes);
}

