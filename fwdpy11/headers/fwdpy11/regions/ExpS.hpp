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
        double mean;

        template <typename Dominance>
        ExpS(const Region& r, double sc, double m, Dominance&& h)
            : Sregion(r, sc, 1, std::forward<Dominance>(h)), mean(m)
        {
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::unique_ptr<ExpS>(
                new ExpS(this->region, this->scaling, this->mean, dominance));
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
                [&rng, this]() {
                    return gsl_ran_exponential(rng.get(), mean) / scaling;
                },
                [this, &rng](const auto esize) {
                    return dominance(rng, esize);
                },
                this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double P) const override
        {
            return gsl_cdf_exponential_Pinv(P, mean) / scaling;
        }
    };
} // namespace fwdpy11

#endif

