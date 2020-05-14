#ifndef FWDPY11_UNIFORMS_HPP
#define FWDPY11_UNIFORMS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct UniformS : public Sregion
    {
        double lo, hi, dominance;

        UniformS(const Region& r, double sc, double lo_, double hi_, double h)
            : Sregion(r, sc, 1), lo(lo_), hi(hi_), dominance(h)
        {
            if (!std::isfinite(lo))
                {
                    throw std::invalid_argument("lo must be finite");
                }
            if (!std::isfinite(hi))
                {
                    throw std::invalid_argument("hi must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
            if (!(hi > lo))
                {
                    throw std::invalid_argument("hi must be > lo");
                }
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::unique_ptr<UniformS>(new UniformS(*this));
        }

        std::uint32_t
        operator()(fwdpp::flagged_mutation_queue& recycling_bin,
                   std::vector<Mutation>& mutations,
                   std::unordered_multimap<double, std::uint32_t>& lookup_table,
                   const std::uint32_t generation, const GSLrng_t& rng) const override
        {
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, false, generation,
                [this, &rng]() { return region(rng); },
                [this, &rng]() { return gsl_ran_flat(rng.get(), lo, hi) / scaling; },
                [this]() { return dominance; }, this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double P) const override
        {
            return gsl_cdf_flat_Pinv(P, lo, hi) / scaling;
        }

        std::vector<double>
        get_dominance() const override
        {
            return {dominance};
        }
    };
} // namespace fwdpy11

#endif

