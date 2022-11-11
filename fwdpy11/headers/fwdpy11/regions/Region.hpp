#ifndef FWDPY11_REGION_HPP
#define FWDPY11_REGION_HPP

#include <gsl/gsl_randist.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <cstdint>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct Region
    {
        double beg, end, weight;
        std::uint16_t label;
        bool coupled;
        Region(double b, double e, double w, bool c, std::uint16_t l)
            : beg(b), end(e), weight((c == true) ? (end - beg) * w : w), label(l),
              coupled(c)
        {
            if (!std::isfinite(beg))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(end))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(weight))
                {
                    throw std::invalid_argument("weight must be finite");
                }
            if (weight < 0.0)
                {
                    throw std::invalid_argument("weight must be non-negative");
                }
            if (!(end > beg))
                {
                    throw std::invalid_argument("end must be greater than beg");
                }
        }

        inline double
        operator()(const GSLrng_t& rng) const
        {
            auto rv = gsl_ran_flat(rng.get(), beg, end);
            while (rv == end)
                {
                    rv = gsl_ran_flat(rng.get(), beg, end);
                }
            return rv;
        }

        inline bool
        valid(double start, double stop) const
        {
            auto beg_ok = beg >= start && beg < stop;
            auto end_ok = end > start && end <= stop;
            return beg_ok && end_ok;
        }
    };
} // namespace fwdpy11

#endif
