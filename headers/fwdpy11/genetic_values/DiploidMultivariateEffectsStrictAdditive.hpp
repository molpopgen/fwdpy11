#ifndef FWDPY11_POP_MULTIVARIATE_STRICT_ADDITIVE_HPP
#define FWDPY11_POP_MULTIVARIATE_STRICT_ADDITIVE_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include "DiploidGeneticValue.hpp"
#include "default_update.hpp"
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>

namespace fwdpy11
{
    struct DiploidMultivariateEffectsStrictAdditive : public DiploidGeneticValue
    {
        std::size_t focal_trait_index;

        DiploidMultivariateEffectsStrictAdditive(std::size_t ndim,
                                                 std::size_t focal_trait,
                                                 const GeneticValueIsTrait *gv2w_,
                                                 const GeneticValueNoise *noise_)
            : DiploidGeneticValue(ndim, gv2w_, noise_), focal_trait_index(focal_trait)
        {
            if (focal_trait_index >= ndim)
                {
                    throw std::invalid_argument(
                        "focal trait index must by < number of traits");
                }
        }

        double
        calculate_gvalue(const fwdpy11::DiploidGeneticValueData data) override
        {
            std::fill(begin(gvalues), end(gvalues), 0.0);

            const auto &pop = data.pop.get();
            const auto diploid_index = data.offspring_metadata.get().label;
            for (auto key :
                 pop.haploid_genomes[pop.diploids[diploid_index].first].smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != gvalues.size())
                        {
                            throw std::runtime_error("dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes), begin(gvalues),
                                   begin(gvalues), std::plus<double>());
                }

            for (auto key :
                 pop.haploid_genomes[pop.diploids[diploid_index].second].smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != gvalues.size())
                        {
                            throw std::runtime_error("dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes), begin(gvalues),
                                   begin(gvalues), std::plus<double>());
                }
            return gvalues[focal_trait_index];
        }

        void
        update(const fwdpy11::DiploidPopulation & /*pop*/) override
        {
        }
    };
} // namespace fwdpy11

#endif
