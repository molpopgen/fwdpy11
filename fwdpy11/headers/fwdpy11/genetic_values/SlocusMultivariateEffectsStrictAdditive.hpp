#ifndef FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP
#define FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include "SlocusPopMultivariateGeneticValueWithMapping.hpp"
#include "default_update.hpp"

namespace fwdpy11
{
    struct SlocusMultivariateEffectsStrictAdditive
        : public SlocusPopMultivariateGeneticValueWithMapping
    {
        std::size_t focal_trait_index;

        SlocusMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const MultivariateGeneticValueToFitnessMap &gv2w_)
            : SlocusPopMultivariateGeneticValueWithMapping(ndim, gv2w_),
              focal_trait_index(focal_trait)
        {
        }

        SlocusMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const MultivariateGeneticValueToFitnessMap &gv2w_,
            const GeneticValueNoise &noise_)
            : SlocusPopMultivariateGeneticValueWithMapping(ndim, gv2w_,
                                                           noise_),
              focal_trait_index(focal_trait)
        {
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const SlocusPop &pop) const
        {
            std::fill(begin(multivariate_effects), end(multivariate_effects),
                      0.0);

            for (auto key :
                 pop.gametes[pop.diploids[diploid_index].first].smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != multivariate_effects.size())
                        {
                            throw std::runtime_error(
                                "dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes),
                                   begin(multivariate_effects),
                                   begin(multivariate_effects),
                                   std::plus<double>());
                }

            for (auto key :
                 pop.gametes[pop.diploids[diploid_index].second].smutations)
                {
                    const auto &mut = pop.mutations[key];
                    if (mut.esizes.size() != multivariate_effects.size())
                        {
                            throw std::runtime_error(
                                "dimensionality mismatch");
                        }
                    std::transform(begin(mut.esizes), end(mut.esizes),
                                   begin(multivariate_effects),
                                   begin(multivariate_effects),
                                   std::plus<double>());
                }
            return multivariate_effects[focal_trait_index];
        }

        pybind11::object
        pickle() const
        {
            return pybind11::make_tuple(multivariate_effects.size(),
                                        focal_trait_index);
        }
    };
} // namespace fwdpy11

#endif
