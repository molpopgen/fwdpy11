#ifndef FWDPY11_POP_MULTIVARIATE_STRICT_ADDITIVE_HPP
#define FWDPY11_POP_MULTIVARIATE_STRICT_ADDITIVE_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include "DiploidPopulationMultivariateGeneticValueWithMapping.hpp"
#include "default_update.hpp"

namespace fwdpy11
{
    struct DiploidMultivariateEffectsStrictAdditive
        : public DiploidPopulationMultivariateGeneticValueWithMapping
    {
        std::size_t focal_trait_index;

        DiploidMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const MultivariateGeneticValueToFitnessMap &gv2w_)
            : DiploidPopulationMultivariateGeneticValueWithMapping(ndim,
                                                                   gv2w_),
              focal_trait_index(focal_trait)
        {
            if (focal_trait_index >= ndim)
                {
                    throw std::invalid_argument(
                        "focal trait index must by < number of traits");
                }
        }

        DiploidMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const MultivariateGeneticValueToFitnessMap &gv2w_,
            const GeneticValueNoise &noise_)
            : DiploidPopulationMultivariateGeneticValueWithMapping(ndim, gv2w_,
                                                                   noise_),
              focal_trait_index(focal_trait)
        {
            if (focal_trait_index >= ndim)
                {
                    throw std::invalid_argument(
                        "focal trait index must by < number of traits");
                }
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const DiploidPopulation &pop) const
        {
            std::fill(begin(gvalues), end(gvalues), 0.0);

            for (auto key :
                 pop.haploid_genomes[pop.diploids[diploid_index].first]
                     .smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != gvalues.size())
                        {
                            throw std::runtime_error(
                                "dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes),
                                   begin(gvalues), begin(gvalues),
                                   std::plus<double>());
                }

            for (auto key :
                 pop.haploid_genomes[pop.diploids[diploid_index].second]
                     .smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != gvalues.size())
                        {
                            throw std::runtime_error(
                                "dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes),
                                   begin(gvalues), begin(gvalues),
                                   std::plus<double>());
                }
            return gvalues[focal_trait_index];
        }

        pybind11::object
        pickle() const
        {
            return pybind11::make_tuple(gvalues.size(), focal_trait_index);
        }
    };
} // namespace fwdpy11

#endif
