#ifndef FWDPY11_SREGION_HPP
#define FWDPY11_SREGION_HPP

#include <memory>
#include "Region.hpp"

namespace fwdpy11
{
    struct Sregion
    {
        Region region; // For returning positions
        std::uint16_t label;
        double scaling;
        Sregion(double b, double e, double w, bool c, std::uint16_t l,
                double s)
            : region(b, e, w, c), label(l), scaling(s)
        {
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

        virtual std::unique_ptr<Sregion> clone() = 0;
    };
} // namespace fwdpy11

#endif
