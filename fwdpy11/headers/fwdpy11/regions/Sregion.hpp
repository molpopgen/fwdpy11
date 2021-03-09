#ifndef FWDPY11_SREGION_HPP
#define FWDPY11_SREGION_HPP

#include <memory>
#include <vector>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/rng.hpp>
#include <gsl/gsl_cdf.h>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>
#include "Region.hpp"

namespace fwdpy11
{
    struct Sregion
    {
        Region region; // For returning positions
        double scaling;
        std::size_t total_dim;
        std::shared_ptr<MutationDominance> dominance;

        virtual ~Sregion() = default;

        template <typename Dominance>
        Sregion(const Region& r, double s, std::size_t dim, Dominance&& h)
            : region(r), scaling(s), total_dim(dim),
              dominance(process_input_dominance(std::forward<Dominance>(h)))
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling must be finite");
                }
            if (dim == 0)
                {
                    throw std::invalid_argument("invalid dimension parameter");
                }
        }

        Sregion(const Sregion&) = delete;
        Sregion(Sregion&&) = default;
        Sregion& operator=(const Sregion&) = delete;
        Sregion& operator=(Sregion&&) = default;

        inline double
        beg() const
        {
            return region.beg;
        }

        inline double
        end() const
        {
            return region.end;
        }

        inline double
        weight() const
        {
            return region.weight;
        }

        inline double
        label() const
        {
            return region.label;
        }

        virtual std::unique_ptr<Sregion> clone() const = 0;
        virtual std::uint32_t
        operator()(fwdpp::flagged_mutation_queue& /*recycling_bin*/,
                   std::vector<Mutation>& /*mutations*/,
                   std::unordered_multimap<double, std::uint32_t>& /*lookup_table*/,
                   const std::uint32_t /*generation*/,
                   const GSLrng_t& /*rng*/) const = 0;
        // Added in 0.7.0.  We now require that these types
        // are able to return deviates from the relevant cdf_P
        // function.
        virtual double from_mvnorm(const double /*deviate*/,
                                   const double /*P*/) const = 0;
        virtual double
        generate_dominance(const GSLrng_t& rng, const double esize) const
        {
            return dominance->generate_dominance(rng, esize);
        }

        virtual std::vector<std::size_t>
        shape() const
        {
            return {total_dim};
        }
    };
} // namespace fwdpy11

#endif
