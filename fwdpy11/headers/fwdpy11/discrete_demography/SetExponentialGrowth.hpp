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
#ifndef FWDPY11_SET_EXPONENTIAL_GROWTH_HPP
#define FWDPY11_SET_EXPONENTIAL_GROWTH_HPP

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <tuple>

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct SetExponentialGrowth
        {
            std::uint32_t when;
            std::int32_t deme;
            double G;
            SetExponentialGrowth(std::uint32_t w, std::int32_t d, double g)
                : when(w), deme(d), G(g)
            {
                if (!std::isfinite(G))
                    {
                        throw std::invalid_argument(
                            "SetExponentialGrowth: growth rate must be "
                            "finite");
                    }
                if (d < 0)
                    {
                        throw std::invalid_argument(
                            "SetExponentialGrowth: deme must be non-negative");
                    }
            }
        };

        inline bool
        operator<(const SetExponentialGrowth& lhs,
                  const SetExponentialGrowth& rhs)
        {
            return std::tie(lhs.when, lhs.deme) < std::tie(rhs.when, rhs.deme);
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
