#ifndef FWDPY11_REGIONS_POISSONINTERVAL_HPP
#define FWDPY11_REGIONS_POISSONINTERVAL_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "GeneticMapUnit.hpp"

namespace fwdpy11
{
    struct PoissonInterval : public GeneticMapUnit
    {
        double beg, end, mean;
        PoissonInterval(double b, double e, double m) : beg(b), end(e), mean(m)
        {
            if (!std::isfinite(b))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(e))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (e <= b)
                {
                    throw std::invalid_argument(
                        "end must be greater than beg");
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
            for (unsigned i = 0; i < n; ++i)
                {
                    breakpoints.push_back(gsl_ran_flat(rng.get(), beg, end));
                }
        }

        pybind11::object
        pickle() const final
        {
            return pybind11::make_tuple(beg, end, mean);
        }

        std::unique_ptr<GeneticMapUnit>
        clone() const final
        {
            return std::unique_ptr<GeneticMapUnit>(new PoissonInterval(*this));
        }

        static PoissonInterval
        unpickle(pybind11::object o)
        {
            auto t = o.cast<pybind11::tuple>();
            if (t.size() != 3)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            return PoissonInterval(t[0].cast<double>(), t[1].cast<double>(),
                                   t[2].cast<double>());
        }
    };
} // namespace fwdpy11

#endif
