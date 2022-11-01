//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWDPY11_PLEIOTROPIC_OPTIMA_HPP
#define FWDPY11_PLEIOTROPIC_OPTIMA_HPP

#include <cmath>
#include <vector>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

namespace fwdpy11
{
    struct PleiotropicOptima
    {
        const std::uint32_t when;
        const std::vector<double> optima;
        const double VW;
        static const std::uint32_t null = std::numeric_limits<std::uint32_t>::max();

        PleiotropicOptima(std::uint32_t w, std::vector<double> o, double vw)
            : when{w}, optima{std::move(o)}, VW{vw}
        {
            if (!std::isfinite(VW))
                {
                    throw std::invalid_argument("VS must be finite");
                }
            if (optima.empty())
                {
                    throw std::invalid_argument("optima cannot be empty");
                }
            if (std::any_of(begin(optima), end(optima),
                            [](double d) { return !std::isfinite(d); })
                == true)
                {
                    throw std::invalid_argument("all optimum values must be finite");
                }
        }

        PleiotropicOptima(std::vector<double> o, double vw)
            : when{std::numeric_limits<std::uint32_t>::max()}, optima{std::move(o)},
              VW{vw}
        {
            if (!std::isfinite(VW))
                {
                    throw std::invalid_argument("VS must be finite");
                }
            if (optima.empty())
                {
                    throw std::invalid_argument("optima cannot be empty");
                }
            if (std::any_of(begin(optima), end(optima),
                            [](double d) { return !std::isfinite(d); })
                == true)
                {
                    throw std::invalid_argument("all optimum values must be finite");
                }
        }
    };

    inline bool
    operator==(const PleiotropicOptima& lhs, const PleiotropicOptima& rhs)
    {
        return lhs.when == rhs.when && lhs.VW == rhs.VW && lhs.optima == rhs.optima;
    }

    inline bool
    operator<(const PleiotropicOptima& lhs, const PleiotropicOptima& rhs)
    {
        return lhs.when < rhs.when;
    }
}

#endif

