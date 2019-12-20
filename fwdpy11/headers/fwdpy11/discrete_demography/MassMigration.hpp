//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef FWDPY11_MASSMIGRATION_HPP
#define FWDPY11_MASSMIGRATION_HPP

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <tuple>

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct MassMigration
        {
            std::int32_t source, destination; // Demes
            std::uint32_t when;               // when it happens
            std::uint32_t number;             // no. to move||copy
            std::int32_t sex;      // If event is sex-specific, this is the sex
            double fraction;       // fraction of deme to move||copy
            bool move_individuals; // If true, we move
            bool
                sex_specific; // If true, only apply to individuals with value sex
            // If size is changed, but not growth rate, do we reset growth rate to zero?
            bool resets_growth_rate;
            MassMigration(std::int32_t s, std::int32_t d, std::uint32_t w,
                          std::uint32_t n, std::int32_t sx, double f, bool mv,
                          bool ss, bool reset)
                : source(s), destination(d), when(w), number(n), sex(sx),
                  fraction(f), move_individuals(mv), sex_specific(ss),
                  resets_growth_rate(reset)
            {
                if (!std::isfinite(fraction))
                    {
                        throw std::invalid_argument(
                            "MassMigration: fraction must be finite");
                    }
                if (fraction <= 0. || fraction > 1.0)
                    {
                        throw std::invalid_argument(
                            "MassMigration: fraction must be 0 < f <= "
                            "1.0");
                    }
                if (source < 0 || destination < 0)
                    {
                        throw std::invalid_argument(
                            "MassMigration: deme indexes must be >=0");
                    }
                if (source == destination)
                    {
                        throw std::invalid_argument(
                            "MassMigration: source must not equal "
                            "destination");
                    }
            }
        };

        inline bool
        operator<(const MassMigration& lhs, const MassMigration& rhs)
        {
            // This operator forces copies to occur before moves
            // from same source to same destination at same time.
            return std::tie(lhs.when, lhs.source, lhs.destination,
                            lhs.move_individuals)
                   < std::tie(rhs.when, rhs.source, rhs.destination,
                              rhs.move_individuals);
        }

        inline bool
        operator==(const MassMigration& lhs, const MassMigration& rhs)
        {
            return lhs.source == rhs.source
                   && lhs.destination == rhs.destination
                   && lhs.when == rhs.when && lhs.number == rhs.number
                   && lhs.fraction == rhs.fraction
                   && lhs.move_individuals == rhs.move_individuals
                   && lhs.resets_growth_rate == rhs.resets_growth_rate;
        }

        inline bool
        operator!=(const MassMigration& lhs, const MassMigration& rhs)
        {
            return !(lhs == rhs);
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
