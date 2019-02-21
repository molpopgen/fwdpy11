#ifndef FWDPY11_REGIONS_BINOMIALPOINT_HPP
#define FWDPY11_REGIONS_BINOMIALPOINT_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "GeneticMapUnit.hpp"

namespace fwdpy11
{
    struct PoissonPoint : public GeneticMapUnit
    {
        double position, mean;
        PoissonPoint(const double pos, const double m) : position(pos), mean(m)
        {
            if (!std::isfinite(pos))
                {
                    throw std::invalid_argument("position must be finite");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (mean < 0)
                {
                    throw std::invalid_argument("mean must be non-negative");
                }
        }

        void
        operator()(const GSLrng_t& rng,
                   std::vector<double>& breakpoints) const final
        {
            unsigned n = gsl_ran_poisson(rng.get(), mean);
            if (n % 2 != 0.0)
                {
                    breakpoints.push_back(position);
                }
        }

        pybind11::object
        pickle() const final
        {
            return pybind11::make_tuple(position, mean);
        }

        static PoissonPoint
        unpickle(pybind11::object o)
        {
            auto t = o.cast<pybind11::tuple>();
            if (t.size() != 2)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            return PoissonPoint(t[0].cast<double>(), t[1].cast<double>());
        }

        std::unique_ptr<GeneticMapUnit>
        clone() const final
        {
            return std::unique_ptr<GeneticMapUnit>(new PoissonPoint(*this));
        }
    };
} // namespace fwdpy11

#endif

