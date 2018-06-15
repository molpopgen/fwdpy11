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
    /// API class
    /// For a single-locus simulation, we need the following concepts:
    /// 1. Calculate the genetic value of a diploid, g.  This is operator()
    /// 2. Calculate any random effects, or "noise", e.  This is nose(...)
    /// 3. Calculate the final fitness, w = f(g,e). This is genetic_value_to_fitness()
    /// 4. If a derived class has any of its own data, provide a means for updating. This is update()
    ///
    /// In a simulation, we want the following operations:
    /// pop.diploid_metadata[i].g = gv(i,pop);
    /// pop.diploid_metadata[i].e = gv.noise(rng, pop.diploid_medatadata[i],p1,p2,pop);
    /// pop.diploid_metadata[i].w = gv.genetic_value_to_fitness(pop.diploid_metadata[i].g,
    ///                                                         pop.diploid_metadata[i].e);
    /// At the end of a generation, update may be called.
    ///
    /// Things to note:
    /// Any deme/geography-specific details must be handled by the derived class.
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
