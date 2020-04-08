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
        const double zeta, sigma, dominance;
        const bool univariate;

        LogNormalS(const Region& r, double scaling, double zeta_, double sigma_,
                   double h)
            : Sregion(r, scaling, 1), zeta(zeta_), sigma(sigma_), dominance(h),
              univariate(true)
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
        LogNormalS(const Region& r, double scaling, double h)
            : Sregion(r, scaling, 1), zeta(std::numeric_limits<double>::quiet_NaN()),
              sigma(std::numeric_limits<double>::quiet_NaN()), dominance(h),
              univariate(false)
        {
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            return std::unique_ptr<LogNormalS>(new LogNormalS(*this));
        }

        std::string
        repr() const override
        {
            std::ostringstream out;
            out.precision(4);
            out << "LogNormalS(";
            this->region.region_repr(out);
            out << ", zeta=" << this->zeta << ", sigma=" << this->sigma
                << ", h=" << this->dominance << ", scaling=" << this->scaling << ')';
            return out.str();
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
                [this]() { return dominance; }, this->label());
        }

        double
        from_mvnorm(const double deviate, const double /*P*/) const override
        {
            return std::exp(deviate) / scaling;
        }

        std::vector<double>
        get_dominance() const override
        {
            return {dominance};
        }

        pybind11::tuple
        pickle() const override
        {
            return pybind11::make_tuple(Sregion::pickle_Sregion(), zeta, sigma,
                                        dominance, univariate);
        }

        static LogNormalS
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 5)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            bool uv = t[4].cast<bool>();
            if (uv)
                {
                    return LogNormalS(Region::unpickle(base[0]), base[1].cast<double>(),
                                      t[1].cast<double>(), t[2].cast<double>(),
                                      t[3].cast<double>());
                }
            return LogNormalS(Region::unpickle(base[0]), base[1].cast<double>(),
                              t[3].cast<double>());
        }
    };
}

#endif

