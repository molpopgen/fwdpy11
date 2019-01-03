#ifndef FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP
#define FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP

#include <vector>
#include "SlocusPopGeneticValueWithMapping.hpp"

namespace fwdpy11
{
    struct SlocusPopMultivariateEffectsStrictAdditive
        : public SlocusPopGeneticValueWithMapping
    {
        mutable std::vector<double> summed_effects;

        SlocusPopMultivariateEffectsStrictAdditive(
            std::size_t ndim, const GeneticValueIsTrait &gv2w_)
            : SlocusPopGeneticValueWithMapping(gv2w_),
              summed_effects(ndim, 0.0)
        {
        }

        SlocusPopMultivariateEffectsStrictAdditive(
            std::size_t ndim, const GeneticValueIsTrait &gv2w_,
            const GeneticValueNoise &noise_)
            : SlocusPopGeneticValueWithMapping(gv2w_, noise_),
              summed_effects(ndim, 0.0)
        {
        }
    };
} // namespace fwdpy11

#endif
