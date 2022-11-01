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
#include <numeric>
#include <limits>
#include <gsl/gsl_matrix.h>
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

            void
            validate_sum(double sum)
            {
                bool close_to_zero = std::fabs(sum - 0.0)
                                     <= 10. * std::numeric_limits<double>::epsilon();
                bool close_to_one = std::fabs(sum - 1.0)
                                    <= 10. * std::numeric_limits<double>::epsilon();
                if (!close_to_one && !close_to_zero)
                    {
                        throw std::invalid_argument(
                            "migration rates must sum to ~0. or ~1. in "
                            "a row.");
                    }
            }

            SetMigrationRates(std::uint32_t w, std::int32_t d, std::vector<double> r)
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

                auto sum = std::accumulate(begin(migrates), end(migrates), 0.0);
                validate_sum(sum);
            }

            SetMigrationRates(std::uint32_t w, std::vector<double> migmatrix)
                : when(w), deme(NULLDEME), migrates(std::move(migmatrix))
            {
                if (migrates.empty())
                    {
                        throw std::invalid_argument("empty migration matrix");
                    }
                double nrow_d = std::sqrt(migrates.size());
                double intpart;
                auto f = std::modf(nrow_d, &intpart);
                if (f != 0.0)
                    {
                        throw std::invalid_argument("input matrix is not square");
                    }
                std::size_t npops{static_cast<std::size_t>(nrow_d)};
                gsl_matrix_const_view v
                    = gsl_matrix_const_view_array(migrates.data(), npops, npops);
                for (std::size_t r = 0; r < npops; ++r)
                    {
                        auto rview = gsl_matrix_const_row(&v.matrix, r);
                        double sum = 0.0;
                        for (std::size_t i = 0; i < rview.vector.size; ++i)
                            {
                                auto v = gsl_vector_get(&rview.vector, i);
                                if (v < 0)
                                    {
                                        throw std::invalid_argument(
                                            "migration rates must be non-negative");
                                    }
                                if (!std::isfinite(v))
                                    {
                                        throw std::invalid_argument(
                                            "migration rates must be finite");
                                    }
                                sum += v;
                            }
                        validate_sum(sum);
                    }
            }
        };

        inline bool
        operator<(const SetMigrationRates& lhs, const SetMigrationRates& rhs)
        {
            return lhs.when < rhs.when;
        }
    } // namespace discrete_demography
} // namespace fwdpy11
#endif
