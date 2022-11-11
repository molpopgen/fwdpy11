#include <cmath>
#include <stdexcept>
#include <functional>
#include <core/internal/gsl_ran_flat.hpp>
#include <gsl/gsl_randist.h> // FIXME: hide in static lib
#include <fwdpp/util/validators.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace fwdpy11
{
    // FIXME: can hide in the static lib
    MutationDominance
    fixed_dominance(double d)
    {
        fwdpp::validators::isfinite(d, "dominance values must be finite");
        return MutationDominance([d](const GSLrng_t&, const double) { return d; });
    }

    MutationDominance
    large_effect_exponentially_recessive(double k, double scaling)
    {
        fwdpp::validators::isfinite(k, "k must be finite");
        fwdpp::validators::isfinite(scaling, "scaling must be finite");
        fwdpp::validators::is_positive(k, "k must be > 0.0");
        return MutationDominance(
            [k, scaling](const GSLrng_t&, const double effect_size) {
                return scaling * std::exp(-k * std::abs(effect_size));
            });
    }

    MutationDominance
    exponential_dominance(double mean)
    {
        fwdpp::validators::isfinite(mean, "mean dominance must be finite");
        return MutationDominance([mean](const GSLrng_t& rng, const double) {
            return gsl_ran_exponential(rng.get(), mean);
        });
    }

    MutationDominance
    uniform_dominance(double lo, double hi)
    {
        fwdpp::validators::isfinite(lo, "lo dominance must be finite");
        fwdpp::validators::isfinite(hi, "hi dominance must be finite");
        if (hi <= lo)
            {
                throw std::invalid_argument("hi must be > lo");
            }
        return MutationDominance([lo, hi](const GSLrng_t& rng, const double) {
            return fwdpy11_core::internal::gsl_ran_flat(rng, lo, hi);
        });
    }

    MutationDominance
    process_input_dominance(double dominance)
    {
        return fixed_dominance(dominance);
    }

    MutationDominance
    process_input_dominance(const MutationDominance& dominance)
    {
        return MutationDominance(dominance);
    }
}
