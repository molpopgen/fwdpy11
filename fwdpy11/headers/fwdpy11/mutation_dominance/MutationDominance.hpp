#ifndef FWDPY11_MUTATION_DOMINANCE_HPP
#define FWDPY11_MUTATION_DOMINANCE_HPP

#include <functional>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct MutationDominance
    {
        using function = std::function<double(const GSLrng_t&, const double)>;
        function generate_dominance;

        MutationDominance(function f) : generate_dominance(std::move(f))
        {
        }

        inline double
        operator()(const GSLrng_t& rng, const double effect_size) const
        {
            return generate_dominance(rng, effect_size);
        }
    };

    MutationDominance fixed_dominance(double d);

    MutationDominance large_effect_exponentially_recessive(double k, double scaling);

    MutationDominance exponential_dominance(double mean);

    MutationDominance uniform_dominance(double lo, double hi);

    MutationDominance process_input_dominance(double dominance);

    MutationDominance process_input_dominance(const MutationDominance& dominance);
}

#endif
