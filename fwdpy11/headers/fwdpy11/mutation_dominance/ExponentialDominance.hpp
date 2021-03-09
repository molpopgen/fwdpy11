#ifndef FWDPY11_EXPONENTIAL_DOMINANCE_HPP
#define FWDPY11_EXPONENTIAL_DOMINANCE_HPP

#include <fwdpp/util/validators.hpp>
#include <gsl/gsl_randist.h>
#include "MutationDominance.hpp"

namespace fwdpy11
{
    struct ExponentialDominance : public MutationDominance
    {
        const double mean;

        explicit ExponentialDominance(double mean) : MutationDominance{}, mean{mean}
        {
            fwdpp::validators::isfinite(mean, "mean dominance must be finite");
        }

        double
        generate_dominance(const GSLrng_t& rng,
                           const double /*effect_size*/) const override final
        {
            return gsl_ran_exponential(rng.get(), mean);
        }

        std::shared_ptr<MutationDominance>
        clone() const override final
        {
            return std::make_shared<ExponentialDominance>(mean);
        }
    };
}

#endif
