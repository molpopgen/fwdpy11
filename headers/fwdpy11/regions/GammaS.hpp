#ifndef FWDPY11_GAMMAS_HPP
#define FWDPY11_GAMMAS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct GammaS : public Sregion
    {
        double mean, shape_parameter;

        template <typename Dominance>
        GammaS(const Region& r, double sc, double m, double s, Dominance&& h)
            : Sregion(r, sc, 1, std::forward<Dominance>(h)), mean(m), shape_parameter(s)
        {
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (!std::isfinite(shape_parameter))
                {
                    throw std::invalid_argument("shape must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::make_unique<GammaS>(this->region, this->scaling, this->mean,
                                            this->shape_parameter, dominance);
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
                    return gsl_ran_gamma(rng.get(), shape_parameter,
                                         mean / shape_parameter)
                           / scaling;
                },
                [this, &rng](const double esize) {
                    return dominance(rng, esize);
                },
                this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double P) const override
        {
            return gsl_cdf_gamma_Pinv(P, shape_parameter, mean / shape_parameter)
                   / scaling;
        }
    };
} // namespace fwdpy11

#endif
