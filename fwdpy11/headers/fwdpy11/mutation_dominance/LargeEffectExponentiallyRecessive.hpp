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

        explicit LargeEffectExponentiallyRecessive(double k) : MutationDominance{}, k{k}
        {
            fwdpp::validators::isfinite(k, "k must be finite");
            fwdpp::validators::is_positive(k, "k must be > 0.0");
        }

        double
        generate_dominance(const GSLrng_t& rng,
                           const double effect_size) const override final
        {
            return std::exp(-k * std::abs(effect_size));
        }

        std::shared_ptr<MutationDominance>
        clone() const override final
        {
            return std::make_shared<LargeEffectExponentiallyRecessive>(*this);
        }
    };
}

#endif

