#ifndef FWDPY11_REGION_HPP
#define FWDPY11_REGION_HPP

#include <cmath>
#include <stdexcept>

namespace fwdpy11
{
    struct Region
    {
        double beg, end, weight;
        bool coupled;
        Region(double b, double e, double w, bool c)
            : beg(b), end(e), weight((c == true) ? (end - beg) * w : w),
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
                    throw std::invalid_argument(
                        "end must be greater than beg");
                }
        }
    };
} // namespace fwdpy11

#endif
