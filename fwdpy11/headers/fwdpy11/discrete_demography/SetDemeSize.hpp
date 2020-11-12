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
#ifndef SET_DEME_SIZE_HPP
#define SET_DEME_SIZE_HPP

#include <cstdint>
#include <stdexcept>

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct SetDemeSize
        {
            std::uint32_t when;
            std::int32_t deme;
            std::uint32_t new_size;
            // If size is changed, but not growth rate, do we reset growth rate to zero?
            bool resets_growth_rate;

            SetDemeSize(std::uint32_t w, std::int32_t d, std::uint32_t n,
                        bool reset)
                : when(w), deme(d), new_size(n), resets_growth_rate(reset)
            {
                if (deme < 0)
                    {
                        throw std::invalid_argument(
                            "SetDemeSize:"
                            " deme must be non-negative");
                    }
            }
        };

        inline bool
        operator<(const SetDemeSize& lhs, const SetDemeSize& rhs)
        {
            return lhs.when < rhs.when;
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif

