#ifndef FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP
#define FWDPY11_SLOCUSPOP_MULTIVARIATE_STRICT_ADDITIVE_HPP

#include <vector>
#include <stdexcept>
#include "SlocusPopGeneticValueWithMapping.hpp"

namespace fwdpy11
{
    struct SlocusMultivariateEffectsStrictAdditive
        : public SlocusPopGeneticValueWithMapping
    {
        mutable std::vector<double> summed_effects;
        std::size_t focal_trait_index;

        SlocusMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const GeneticValueIsTrait &gv2w_)
            : SlocusPopGeneticValueWithMapping(gv2w_),
              summed_effects(ndim, 0.0), focal_trait_index(focal_trait)
        {
        }

        SlocusMultivariateEffectsStrictAdditive(
            std::size_t ndim, std::size_t focal_trait,
            const GeneticValueIsTrait &gv2w_, const GeneticValueNoise &noise_)
            : SlocusPopGeneticValueWithMapping(gv2w_, noise_),
              summed_effects(ndim, 0.0), focal_trait_index(focal_trait)
        {
        }

        void
        operator()(const GSLrng_t &rng, std::size_t diploid_index,
                   const SlocusPop &pop, DiploidMetadata &metadata) const
        {
            metadata.g = calculate_gvalue(diploid_index, pop);
            metadata.e = noise_fxn->operator()(
                rng, metadata, metadata.parents[0], metadata.parents[1], pop);
            metadata.w = gv2w->operator()(metadata);
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const SlocusPop &pop) const
        {
            //TODO: implement!
            return 0.0;
        }

        void
        update(const SlocusPop &pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }

        pybind11::object
        pickle() const
        {
            return pybind11::none();
        }
    };
} // namespace fwdpy11

#endif
