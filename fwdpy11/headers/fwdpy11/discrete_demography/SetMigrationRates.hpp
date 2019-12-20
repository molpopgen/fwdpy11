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

#ifndef FWDPY11_SET_MIGRATION_RATES_HPP
#define FWDPY11_SET_MIGRATION_RATES_HPP

#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>
#include <limits>
#include <pybind11/numpy.h>
#include "constants.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct SetMigrationRates
        {
            std::uint32_t when;
            std::int32_t deme; // source deme
            std::vector<double> migrates;
            SetMigrationRates(std::uint32_t w, std::int32_t d,
                              std::vector<double> r)
                : when(w), deme(d), migrates(std::move(r))
            {
                if (deme < 0)
                    {
                        throw std::invalid_argument(
                            "SetMigrationRates: deme label must be "
                            "non-negative");
                    }
                if (migrates.empty())
                    {
                        throw std::invalid_argument(
                            "SetMigrationRates: empty list of rates");
                    }
                for (auto r : migrates)
                    {
                        if (r < 0.0)
                            {
                                throw std::invalid_argument(
                                    "SetMigrationRates: rates must be "
                                    "non-negative");
                            }
                        if (!std::isfinite(r))
                            {
                                throw std::invalid_argument(
                                    "SetMigrationRates: rates must be finite");
                            }
                    }
            }
            SetMigrationRates(std::uint32_t w, pybind11::array_t<double> m)
                : when(w), deme(NULLDEME), migrates()
            {
                auto r = m.unchecked<2>();
                if (r.shape(0) != r.shape(1))
                    {
                        throw std::invalid_argument(
                            "SetMigrationRates: input matrix must be square");
                    }
                for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                    {
                        for (decltype(i) j = 0; j < r.shape(1); ++j)
                            {
                                if (r(i, j) < 0.0)
                                    {
                                        throw std::invalid_argument(
                                            "SetMigrationRates: rates must be "
                                            "non-negative");
                                    }
                                if (!std::isfinite(r(i, j)))
                                    {
                                        throw std::invalid_argument(
                                            "SetMigrationRates: rates must be "
                                            "finite");
                                    }
                                migrates.push_back(r(i, j));
                            }
                    }
            }
        };

        inline bool
        operator<(const SetMigrationRates& lhs, const SetMigrationRates& rhs)
        {
            return std::tie(lhs.when, lhs.deme) < std::tie(rhs.when, rhs.deme);
        }
    } // namespace discrete_demography
} // namespace fwdpy11
#endif
