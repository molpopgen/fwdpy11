#pragma once

#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include "fwdpp_wrappers/fwdpp_genetic_value.hpp"

namespace fwdpy11
{
    struct single_deme_additive_het
    {
        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += m.s * m.h;
        }
    };

    struct multi_deme_additive_het
    {
        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += m.esizes[deme] * m.heffects[deme];
        }
    };

    struct single_deme_additive_hom
    {
        double scaling;
        single_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.s;
        }
    };

    struct multi_deme_additive_hom
    {
        double scaling;
        multi_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.esizes[deme];
        }
    };

    struct final_additive_trait
    {
        inline double
        operator()(double d) const
        {
            return d;
        }
    };

    struct final_additive_fitness
    {
        inline double
        operator()(double d) const
        {
            return std::max(0.0, 1. + d);
        }
    };

    using DiploidAdditive = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
        single_deme_additive_het, single_deme_additive_hom, multi_deme_additive_het,
        multi_deme_additive_hom, 0>;

    inline DiploidAdditive
    additive_fitness_model(std::size_t ndemes, double scaling,
                           const GeneticValueNoise* noise)
    {
        return DiploidAdditive(
            ndemes, scaling, final_additive_fitness(),
            [](const double) { return false; }, nullptr, noise);
    }

    inline DiploidAdditive
    additive_trait_model(std::size_t ndemes, double scaling,
                         const GeneticValueIsTrait* gvalue_to_fitness,
                         const GeneticValueNoise* noise)
    {
        return DiploidAdditive(
            ndemes, scaling, final_additive_trait(), [](const double) { return false; },
            gvalue_to_fitness, noise);
    }

}
