#pragma once

#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include "fwdpp_wrappers/fwdpp_genetic_value.hpp"

namespace fwdpy11
{
    struct single_deme_multiplicative_het
    {
        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + m.s * m.h);
        }
    };

    struct multi_deme_multiplicative_het
    {
        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + m.esizes[deme] * m.heffects[deme]);
        }
    };

    struct single_deme_multiplicative_hom
    {
        double scaling;
        single_deme_multiplicative_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + scaling * m.s);
        }
    };

    struct multi_deme_multiplicative_hom
    {
        double scaling;
        multi_deme_multiplicative_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + scaling * m.esizes[deme]);
        }
    };

    struct final_multiplicative_trait
    {
        inline double
        operator()(double d) const
        {
            return d - 1.0;
        }
    };

    struct final_multiplicative_fitness
    {
        inline double
        operator()(double d) const
        {
            return std::max(0.0, d);
        }
    };

    using DiploidMultiplicative
        = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
            single_deme_multiplicative_het, single_deme_multiplicative_hom,
            multi_deme_multiplicative_het, multi_deme_multiplicative_hom, 1>;

    inline DiploidMultiplicative
    multiplicative_fitness_model(std::size_t ndemes, double scaling,
                                 const GeneticValueNoise* noise)
    {
        return DiploidMultiplicative(
            ndemes, scaling, final_multiplicative_fitness(),
            [](const double w) { return w <= 0.0; }, nullptr, noise);
    }

    inline DiploidMultiplicative
    multiplicative_trait_model(std::size_t ndemes, double scaling,
                               const GeneticValueIsTrait* gvalue_to_fitness,
                               const GeneticValueNoise* noise)
    {
        return DiploidMultiplicative(
            ndemes, scaling, final_multiplicative_trait(),
            [](const double) { return false; }, gvalue_to_fitness, noise);
    }
}
