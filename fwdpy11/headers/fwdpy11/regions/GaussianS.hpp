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
        double sd, dominance;
        GaussianS(double b, double e, double w, double sd, double h, bool c,
                  std::uint16_t l, double sc)
            : Sregion(b, e, w, c, l, sc), sd(sd), dominance(h)
        {
            if (!std::isfinite(sd))
                {
                    throw std::invalid_argument("sd must be finite");
                }
            if (!(sd > 0))
                {
                    throw std::invalid_argument("sd must be > 0");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<GaussianS>(new GaussianS(*this));
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
                    return gsl_ran_gaussian_ziggurat(rng.get(), sd) / scaling;
                },
                [this]() { return dominance; }, this->label());
        }
    };
} // namespace fwdpy11

#endif

