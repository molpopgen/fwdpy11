#ifndef FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__
#define FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__

#include <memory>
#include <cstdint>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include "GeneticValueToFitness.hpp"
#include "noise.hpp"

namespace fwdpy11
{
    struct SlocusPopGeneticValue
    ///API class
    {
        virtual double operator()(const std::size_t /*diploid_index*/,
                                  const SlocusPop& /*pop*/) const = 0;
        virtual double genetic_value_to_fitness(const double /*g*/,
                                                const double /*e*/) const = 0;
        virtual double noise(const GSLrng_t& /*rng*/,
                             const dip_metadata& /*offspring_metadata*/,
                             const std::size_t /*parent1*/,
                             const std::size_t /*parent2*/,
                             const SlocusPop& /*pop*/) const = 0;
        virtual void update(const SlocusPop& /*pop*/) = 0;
    };

    struct SlocusPopGeneticValueWithMapping : public SlocusPopGeneticValue
    ///API class.
    {
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<GeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<SlocusPopGeneticValueNoise> noise_fxn;

        SlocusPopGeneticValueWithMapping(
            std::unique_ptr<GeneticValueToFitnessMap> gv2w_)
            : gv2w{ std::move(gv2w_) }, noise_fxn{ new SlocusPopNoNoise() }
        {
        }

        SlocusPopGeneticValueWithMapping(
            std::unique_ptr<GeneticValueToFitnessMap> gv2w_,
            std::unique_ptr<SlocusPopGeneticValueNoise> noise_)
            : gv2w{ std::move(gv2w_) }, noise_fxn{ std::move(noise_) }
        {
        }

        inline virtual double
        genetic_value_to_fitness(const double g, const double e) const
        {
            return gv2w->operator()(g, e);
        }

        inline virtual double
        noise(const GSLrng_t& rng, const dip_metadata& offspring_metadata,
              const std::size_t parent1, const std::size_t parent2,
              const SlocusPop& pop) const
        {
            return noise_fxn->operator()(rng, offspring_metadata, parent1,
                                         parent2, pop);
        }
    };
} //namespace fwdpy11

#endif
