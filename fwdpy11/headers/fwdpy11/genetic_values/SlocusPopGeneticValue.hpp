#ifndef FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__
#define FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__

#include <memory>
#include <cstdint>
#include <fwdpy11/types/SlocusPop.hpp>
#include "GeneticValueToFitness.hpp"

namespace fwdpy11
{
    struct SlocusPopGeneticValue
    ///API class
    {
        virtual double operator()(const std::size_t /*diploid_index*/,
                                         const SlocusPop& /*pop*/) const = 0;
        virtual double
        genetic_value_to_fitness(const double /*g*/,
                                 const double /*e*/) const = 0;
        virtual void update(const SlocusPop& /*pop*/) = 0;
    };

    struct SlocusPopGeneticValueWithMapping : public SlocusPopGeneticValue
    ///API class.
    {
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<GeneticValueToFitnessMap> gv2w;
        SlocusPopGeneticValueWithMapping(
            std::unique_ptr<GeneticValueToFitnessMap> gv2w_)
            : gv2w{ std::move(gv2w_) }
        {
        }

        inline virtual double
        genetic_value_to_fitness(const double g, const double e) const
        {
            return gv2w->operator()(g, e);
        }
    };
} //namespace fwdpy11

#endif
