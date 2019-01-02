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
        UniformS(double b, double e, double w, double lo_, double hi_,
                 double h, bool c, std::uint16_t l, double sc)
            : Sregion(b, e, w, c, l, sc), lo(lo_), hi(hi_), dominance(h)
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
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<UniformS>(new UniformS(*this));
        }

        std::uint32_t
        operator()(
            fwdpp::flagged_mutation_queue& recycling_bin,
            std::vector<Mutation>& mutations,
            std::unordered_multimap<double, std::uint32_t>& lookup_table,
            const std::uint32_t generation, const GSLrng_t& rng) const
        {
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, generation,
                [this, &rng]() { return region(rng); },
                [this, &rng]() { return gsl_ran_flat(rng.get(), lo, hi); },
                [this]() { return dominance; }, this->label());
        }
    };
} // namespace fwdpy11

#endif

