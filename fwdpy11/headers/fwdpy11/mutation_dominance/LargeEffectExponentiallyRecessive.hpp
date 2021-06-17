#ifndef FWDPY11_LARGEEFFECTEXPONENTIALLYRECESSIVE_HPP
#define FWDPY11_LARGEEFFECTEXPONENTIALLYRECESSIVE_HPP

#include <cmath>
#include <fwdpp/util/validators.hpp>
#include "MutationDominance.hpp"

namespace fwdpy11
{
    struct LargeEffectExponentiallyRecessive : public MutationDominance
    {
        const double k;
        const double scaling;

        explicit LargeEffectExponentiallyRecessive(double k, double scaling)
            : MutationDominance{}, k{k}, scaling{scaling}
        {
            fwdpp::validators::isfinite(k, "k must be finite");
            fwdpp::validators::isfinite(scaling, "scaling must be finite");
            fwdpp::validators::is_positive(k, "k must be > 0.0");
        }

        double
        generate_dominance(const GSLrng_t& /*rng*/,
                           const double effect_size) const override final
        {
            return scaling * std::exp(-k * std::abs(effect_size));
        }

        std::shared_ptr<MutationDominance>
        clone() const override final
        {
            return std::make_shared<LargeEffectExponentiallyRecessive>(k, scaling);
        }
    };
}

#endif

