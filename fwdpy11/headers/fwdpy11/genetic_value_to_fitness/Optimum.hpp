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
#ifndef FWDPY11_OPTUMUM_HPP
#define FWDPY11_OPTUMUM_HPP

#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace fwdpy11
{
    struct Optimum
    {
        const std::uint32_t when;
        const double opt;
        const double
            VW; // NOTE: this is used as VS until we sort out nomenclature issues
        static const std::uint32_t null = std::numeric_limits<std::uint32_t>::max();

        Optimum(std::uint32_t w, double o, double vw) : when(w), opt(o), VW(vw)
        {
            if (!std::isfinite(opt) || !std::isfinite(VW))
                {
                    throw std::invalid_argument("opt and VS must both be finite");
                }
            if (VW <= 0.0)
                {
                    throw std::invalid_argument("VS must be >= 0.0");
                }
        }

        Optimum(double o, double vw) : when(null), opt(o), VW(vw)
        {
            if (!std::isfinite(opt) || !std::isfinite(VW))
                {
                    throw std::invalid_argument("opt and VS must both be finite");
                }
            if (VW <= 0.0)
                {
                    throw std::invalid_argument("VS must be >= 0.0");
                }
        }
    };

}

#endif
