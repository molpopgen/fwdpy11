#ifndef FWDPY11_REGIONS_BINOMIALPOINT_HPP
#define FWDPY11_REGIONS_BINOMIALPOINT_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "GeneticMapUnit.hpp"

namespace fwdpy11
{
    struct BinomialPoint : public GeneticMapUnit
    {
        double position, prob;
        BinomialPoint(const double pos, const double pr)
            : position(pos), prob(pr)
        {
            if (!std::isfinite(pos))
                {
                    throw std::invalid_argument("position must be finite");
                }
            if (pr < 0 || pr > 1.)
                {
                    throw std::invalid_argument(
                        "probability must be 0 <= x <= 1");
                }
        }

        void
        operator()(const GSLrng_t& rng,
                   std::vector<double>& breakpoints) const final
        {
            if (gsl_rng_uniform(rng.get()) <= prob)
                {
                    breakpoints.push_back(position);
                }
        }

        pybind11::object
        pickle() const final
        {
            return pybind11::make_tuple(position, prob);
        }

        static BinomialPoint
        unpickle(pybind11::object o)
        {
            auto t = o.cast<pybind11::tuple>();
            if (t.size() != 2)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            return BinomialPoint(t[0].cast<double>(), t[1].cast<double>());
        }

        std::unique_ptr<GeneticMapUnit>
        clone() const final
        {
            return std::unique_ptr<GeneticMapUnit>(new BinomialPoint(*this));
        }
    };
} // namespace fwdpy11

#endif
