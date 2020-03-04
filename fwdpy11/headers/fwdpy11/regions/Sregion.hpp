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
#include "Region.hpp"

namespace fwdpy11
{
    struct Sregion
    {
        Region region; // For returning positions
        double scaling;
        const std::size_t total_dim;

        virtual ~Sregion() = default;

        Sregion(const Region& r, double s, std::size_t dim)
            : region(r), scaling(s), total_dim(dim)
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
        virtual std::string repr() const = 0;
        virtual pybind11::tuple pickle() const = 0;
        virtual std::uint32_t operator()(
            fwdpp::flagged_mutation_queue& /*recycling_bin*/,
            std::vector<Mutation>& /*mutations*/,
            std::unordered_multimap<double, std::uint32_t>& /*lookup_table*/,
            const std::uint32_t /*generation*/,
            const GSLrng_t& /*rng*/) const = 0;
        // Added in 0.7.0.  We now require that these types
        // are able to return deviates from the relevant cdf_P
        // function.
        virtual double from_mvnorm(const double /*deviate*/,
                                   const double /*P*/) const = 0;
        virtual std::vector<double> get_dominance() const = 0;
        virtual pybind11::tuple
        shape() const
        {
            return pybind11::make_tuple(total_dim);
        }

        pybind11::tuple
        pickle_Sregion() const
        {
            return pybind11::make_tuple(region.pickle(), scaling);
        }

        static std::tuple<Region, double>
        unpickle_Sregion(pybind11::tuple t)
        {
            if (t.size() != 2)
                {
                    throw std::runtime_error("inalid tuple size");
                }
            return std::make_tuple(Region::unpickle(t[0]),
                                   t[1].cast<double>());
        }

        inline bool
        is_equal(const Sregion& rhs) const
        {
            return this->region == rhs.region && this->scaling == rhs.scaling;
        }
    };
} // namespace fwdpy11

#endif
