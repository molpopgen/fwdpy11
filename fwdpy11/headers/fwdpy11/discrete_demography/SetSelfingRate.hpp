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
#ifndef FWDPY11_SET_SELFING_RATE_HPP
#define FWDPY11_SET_SELFING_RATE_HPP

#include <cstdint>
#include <cmath>
#include <stdexcept>

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct SetSelfingRate
        {
            std::uint32_t when;
            std::int32_t deme;
            double S;
            SetSelfingRate(std::uint32_t w, std::int32_t d, double s)
                : when(w), deme(d), S(s)
            {
                if (!std::isfinite(S))
                    {
                        throw std::invalid_argument(
                            "SetSelfingRate: selfing rate must be finite");
                    }
                if (S < 0. || S > 1.)
                    {
                        throw std::invalid_argument(
                            "SetSelfingRate: selfing rate must be 0 <= S <= "
                            "1.");
                    }
                if (d < 0)
                    {
                        throw std::invalid_argument(
                            "SetSelfingRate: deme must be non-negative");
                    }
            }
        };

        inline bool
        operator<(const SetSelfingRate& lhs, const SetSelfingRate& rhs)
        {
            return lhs.when < rhs.when;
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif

