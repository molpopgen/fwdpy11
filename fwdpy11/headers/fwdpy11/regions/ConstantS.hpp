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
        bool is_neutral;

        ConstantS(const Region& r, const double s, const double es, const double h)
            : Sregion(r, s, 1), esize(es), dominance(h), is_neutral(false)
        {
            if (!std::isfinite(esize))
                {
                    throw std::invalid_argument("esize must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
            if (esize == 0.0)
                {
                    throw std::invalid_argument("effect size cannot be 0.0");
                }
        }

        // Constructor added in 0.6.3 to allow the back-end
        // of MutationRegions to supply neutral variants.
        // This constructor is NOT exposed to Python.
        ConstantS(const Region& r)
            : Sregion(r, 1., 1), esize(0.), dominance(1.), is_neutral(true)
        {
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::unique_ptr<ConstantS>(new ConstantS(*this));
        }

        std::uint32_t
        operator()(fwdpp::flagged_mutation_queue& recycling_bin,
                   std::vector<Mutation>& mutations,
                   std::unordered_multimap<double, std::uint32_t>& lookup_table,
                   const std::uint32_t generation, const GSLrng_t& rng) const override
        {
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, is_neutral, generation,
                [this, &rng]() { return region(rng); },
                [this]() { return esize / scaling; }, [this]() { return dominance; },
                this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double /*P*/) const override
        {
            //NOTE: ignores the input!
            return esize / scaling;
        }

        std::vector<double>
        get_dominance() const override
        {
            return {dominance};
        }
    };
} // namespace fwdpy11

#endif

