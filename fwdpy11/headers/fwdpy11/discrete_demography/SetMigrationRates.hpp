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

#include <gsl/gsl_matrix.h>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <vector>
#include <numeric>
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
                auto sum
                    = std::accumulate(begin(migrates), end(migrates), 0.0);
                if (sum != 0. && sum != 1.0)
                    {
                        throw std::invalid_argument(
                            "migration rates must sum to 0. or 1. in "
                            "a row.");
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
                std::size_t npops = r.shape(0);
                gsl_matrix_const_view v = gsl_matrix_const_view_array(
                    migrates.data(), npops, npops);
                for (std::size_t i = 0; i < npops; ++i)
                    {
                        gsl_vector_const_view r
                            = gsl_matrix_const_row(&v.matrix, i);
                        double sum = 0.0;
                        for (std::size_t j = 0; j < npops; ++j)
                            {
                                sum += gsl_vector_get(&r.vector, j);
                            }
                        if (sum != 0. && sum != 1.0)
                            {
                                throw std::invalid_argument(
                                    "migration rates must sum to 0. or 1. in "
                                    "a row.");
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
