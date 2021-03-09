#ifndef FWDPY11_REGIONS_LOGNORMALS_HPP__
#define FWDPY11_REGIONS_LOGNORMALS_HPP__

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{
    struct LogNormalS : public Sregion
    {
        const double zeta, sigma;
        const bool univariate;

        template <typename Dominance>
        LogNormalS(const Region& r, double scaling, double zeta_, double sigma_,
                   Dominance&& h)
            : Sregion(r, scaling, 1, std::forward<Dominance>(h)), zeta(zeta_),
              sigma(sigma_), univariate(true)
        {
            if (!std::isfinite(zeta))
                {
                    throw std::invalid_argument("zeta must be finite");
                }
            if (!std::isfinite(sigma))
                {
                    throw std::invalid_argument("sigma must be finite");
                }
            if (sigma <= 0)
                {
                    throw std::invalid_argument("sigma must be > 0");
                }
        }

        // For use with mvS
        template <typename Dominance>
        LogNormalS(const Region& r, double scaling, Dominance&& h)
            : Sregion(r, scaling, 1, std::forward<Dominance>(h)),
              zeta(std::numeric_limits<double>::quiet_NaN()),
              sigma(std::numeric_limits<double>::quiet_NaN()), univariate(false)
        {
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            if (std::isfinite(zeta))
                {
                    return std::unique_ptr<LogNormalS>(
                        new LogNormalS(region, scaling, zeta, sigma, dominance));
                }
            return std::unique_ptr<LogNormalS>(
                new LogNormalS(region, scaling, dominance));
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
                    return gsl_ran_lognormal(rng.get(), zeta, sigma) / scaling;
                },
                [this, &rng](const double esize) {
                    return dominance->generate_dominance(rng, esize);
                },
                this->label());
        }

        double
        from_mvnorm(const double deviate, const double /*P*/) const override
        {
            return std::exp(deviate) / scaling;
        }
    };
}

#endif

