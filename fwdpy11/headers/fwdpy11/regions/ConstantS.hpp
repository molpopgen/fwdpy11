#ifndef FWDPY11_CONSTANTS_HPP
#define FWDPY11_CONSTANTS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct ConstantS : public Sregion
    {
        double esize, dominance;
        ConstantS(double b, double e, double w, double es, double h, bool c,
                  std::uint16_t l, double sc)
            : Sregion(b, e, w, c, l, sc), esize(es), dominance(h)
        {
            if (!std::isfinite(esize))
                {
                    throw std::invalid_argument("esize must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<ConstantS>(new ConstantS(*this));
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
                [this]() { return esize/scaling; },
                [this]() { return dominance; },
                this->label());
        }
    };
} // namespace fwdpy11

#endif

