#ifndef FWDPY11_UNIFORM_DOMINANCE_HPP
#define FWDPY11_UNIFORM_DOMINANCE_HPP

#include <fwdpp/util/validators.hpp>
#include <gsl/gsl_randist.h>
#include <stdexcept>
#include "MutationDominance.hpp"

namespace fwdpy11
{
    struct UniformDominance : public MutationDominance
    {
        const double lo, hi;

        UniformDominance(double lo, double hi) : MutationDominance{}, lo{lo}, hi{hi}
        {
            fwdpp::validators::isfinite(lo, "lo dominance must be finite");
            fwdpp::validators::isfinite(hi, "hi dominance must be finite");
            if (hi <= lo)
                {
                    throw std::invalid_argument("hi must be > lo");
                }
        }

        double
        generate_dominance(const GSLrng_t& rng,
                           const double /*effect_size*/) const override final
        {
            return gsl_ran_flat(rng.get(), lo, hi);
        }

        std::shared_ptr<MutationDominance>
        clone() const override final
        {
            return std::make_shared<UniformDominance>(lo, hi);
        }
    };
}

#endif

