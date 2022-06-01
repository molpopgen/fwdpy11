#ifndef FWDPY11_GAUSSIANS_HPP
#define FWDPY11_GAUSSIANS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct GaussianS : public Sregion
    {
        double sd;

        template <typename Dominance>
        GaussianS(const Region& r, double sc, double sd_, Dominance&& h)
            : Sregion(r, sc, 1, std::forward<Dominance>(h)), sd(sd_)
        {
            if (!std::isfinite(sd))
                {
                    throw std::invalid_argument("sd must be finite");
                }
            if (!(sd > 0))
                {
                    throw std::invalid_argument("sd must be > 0");
                }
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::make_unique<GaussianS>(this->region, this->scaling, this->sd,
                                               dominance);
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
                [this, &rng]() {
                    return gsl_ran_gaussian_ziggurat(rng.get(), sd) / scaling;
                },
                [this, &rng](const double esize) {
                    return dominance(rng, esize);
                },
                this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double P) const override
        {
            return gsl_cdf_gaussian_Pinv(P, sd) / scaling;
        }
    };
} // namespace fwdpy11

#endif

