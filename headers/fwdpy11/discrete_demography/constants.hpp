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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_CONSTANTS_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_CONSTANTS_HPP

#include <cstdint>
#include <limits>

namespace fwdpy11
{
    namespace discrete_demography
    {
        constexpr double NOGROWTH = 1.0;
        constexpr std::int32_t NULLDEME
            = std::numeric_limits<std::int32_t>::min();
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
