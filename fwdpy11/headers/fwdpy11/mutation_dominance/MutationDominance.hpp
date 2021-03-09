#ifndef FWDPY11_MUTATION_DOMINANCE_HPP
#define FWDPY11_MUTATION_DOMINANCE_HPP

#include <cmath>
#include <memory>
#include <stdexcept>
#include <fwdpp/util/validators.hpp>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct MutationDominance
    {
        virtual ~MutationDominance() = default;
        MutationDominance() = default;
        MutationDominance(const MutationDominance&) = delete;
        MutationDominance(MutationDominance&&) = default;
        MutationDominance& operator=(const MutationDominance&) = delete;
        MutationDominance& operator=(MutationDominance&&) = default;
        virtual double generate_dominance(const GSLrng_t& /*rng*/,
                                          const double /*effect_size*/) const = 0;
        virtual std::shared_ptr<MutationDominance> clone() const = 0;
    };

    struct FixedDominance : public MutationDominance
    {
        double dominance;
        explicit FixedDominance(double d) : dominance{d}
        {
            fwdpp::validators::isfinite(d, "dominance values must be finite");
        }

        double
        generate_dominance(const GSLrng_t& /*rng*/,
                           const double /*effect_size*/) const override final
        {
            return dominance;
        }

        std::shared_ptr<MutationDominance>
        clone() const override final
        {
            return std::make_shared<FixedDominance>(dominance);
        }
    };

    inline std::shared_ptr<MutationDominance>
    process_input_dominance(double dominance)
    {
        return std::make_shared<FixedDominance>(dominance);
    }

    // NOTE/FIXME we may not need both of these overloads:

    inline std::shared_ptr<MutationDominance>
    process_input_dominance(const MutationDominance& dominance)
    {
        return dominance.clone();
    }

    inline std::shared_ptr<MutationDominance>
    process_input_dominance(const std::shared_ptr<MutationDominance>& dominance)
    {
        return dominance->clone();
    }
}

#endif
