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
#ifndef FWDPY11_GENETIC_VALUES_DETAILS_MLOCUS_AGGREGATORS_HPP__
#define FWDPY11_GENETIC_VALUES_DETAILS_MLOCUS_AGGREGATORS_HPP__

#include <vector>
#include <numeric>
#include <functional>

namespace fwdpy11
{
    struct aggregate_additive_fitness
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::max(0., std::accumulate(g.data(), g.data() + s, 0.0)
                                    - (s - 1));
        }
    };

    struct aggregate_mult_fitness
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::max(0., std::accumulate(g.data(), g.data() + s, 1.0,
                                                std::multiplies<double>()));
        }
    };

    struct aggregate_additive_trait
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            return std::accumulate(g.data(), g.data() + g.size(), 0.0);
        }
    };

    struct aggregate_mult_trait
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::accumulate(
                       g.data(), g.data() + s, 1.0,
                       [](double prod, double v) { return prod * (1. + v); })
                   - 1.0;
        }
    };
} // namespace fwdpy11

#endif
