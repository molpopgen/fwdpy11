#ifndef FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__
#define FWDPY11_SLOCUSPOP_GENETIC_VALUE_HPP__

#include <cstdint>
#include <fwdpy11/types/SlocusPop.hpp>

namespace fwdpy11
{
    struct SlocusPopGeneticValue
    ///API class
    {
        inline virtual double operator()(const std::size_t /*diploid_index*/,
                                         const SlocusPop& /*pop*/) const = 0;
        inline virtual double
        genetic_value_to_fitness(const double /*g*/,
                                 const double /*e*/) const = 0;
        inline virtual void update(const SlocusPop& /*pop*/) = 0;
    };
} //namespace fwdpy11

#endif
