#ifndef FWDPY11_EXPS_HPP
#define FWDPY11_EXPS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct ExpS : public Sregion
    {
        double mean, dominance;
        ExpS(double b, double e, double w, double m, double h, bool c,
             std::uint16_t l, double sc)
            : Sregion(b, e, w, c, l, sc), mean(m), dominance(h)
        {
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<ExpS>(new ExpS(*this));
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
                [this, &rng]() {
                    return gsl_ran_exponential(rng.get(), mean);
                },
                [this]() { return dominance; }, this->label);
        }
    };
} // namespace fwdpy11

#endif

