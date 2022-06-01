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
        double esize;
        bool is_neutral;

        template <typename Dominance>
        ConstantS(const Region& r, const double s, const double es, Dominance&& h)
            : Sregion(r, s, 1, std::forward<Dominance>(h)), esize(es), is_neutral(false)
        {
            if (!std::isfinite(esize))
                {
                    throw std::invalid_argument("esize must be finite");
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
            : Sregion(r, 1., 1, process_input_dominance(1.)), esize(0.), is_neutral(true)
        {
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::make_unique<ConstantS>(this->region, this->scaling, this->esize,
                                               dominance);
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
                [this]() { return esize / scaling; },
                [this, &rng](const double esize) { return dominance(rng, esize); },
                this->label());
        }

        double
        from_mvnorm(const double /*deviate*/, const double /*P*/) const override
        {
            //NOTE: ignores the input!
            return esize / scaling;
        }
    };
} // namespace fwdpy11

#endif
