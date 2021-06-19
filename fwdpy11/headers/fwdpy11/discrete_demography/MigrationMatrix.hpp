
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
#ifndef FWDPY11_MIGRATION_MATRIX_HPP
#define FWDPY11_MIGRATION_MATRIX_HPP

#include <vector>
#include <memory>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <fwdpp/gsl_discrete.hpp>
#include "constants.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class MigrationMatrix
        {
          private:
            void
            validate_row(std::size_t i) const
            {
                if (i * npops >= M.size())
                    {
                        throw std::invalid_argument(
                            "MigrationMatrix: row index out of range");
                    }
                double rsum = 0.;
                for (std::size_t j = i * npops; j < i * npops + npops; ++j)
                    {
                        auto v = M[j];
                        if (v < 0.0)
                            {
                                throw std::invalid_argument(
                                    "migration rates must be non-negative");
                            }
                        if (!std::isfinite(v))
                            {
                                throw std::invalid_argument(
                                    "migration rates must be finite");
                            }
                        rsum += v;
                    }
                if (rsum != 0. && rsum != 1.0)
                    {
                        throw std::invalid_argument(
                            "migration rates must sum to 0. or 1. in a row.");
                    }
            }

            void
            validate_all_row_sums() const
            {
                for (std::size_t i = 0; i < npops; ++i)
                    {
                        validate_row(i);
                    }
            }

          public:
            std::vector<double> M;
            std::size_t npops;
            bool scaled;

            MigrationMatrix(): M{}, npops{0}, scaled{false} {}

            template <typename T>
            MigrationMatrix(T&& matrix, std::size_t nrows,
                            const bool scaled_rates)
                : M(std::forward<T>(matrix)), npops(nrows),
                  scaled(scaled_rates)
            {
                if (M.size() % nrows != 0.0)
                    {
                        throw std::invalid_argument(
                            "MigrationMatrix must be square");
                    }
                validate_all_row_sums();
            }

            std::unique_ptr<MigrationMatrix>
            clone() const
            {
                return std::unique_ptr<MigrationMatrix>(
                    new MigrationMatrix(*this));
            }

            void
            set_migration_rates(std::int32_t source,
                                const std::vector<double>& migrates)
            {
                if (source != NULLDEME
                    && static_cast<std::size_t>(source) >= npops)
                    {
                        throw std::invalid_argument(
                            "source pop index out of range");
                    }
                if (source != NULLDEME && migrates.size() != npops)
                    {
                        throw std::invalid_argument(
                            "invalid number of migration rates");
                    }
                else if (source == NULLDEME && migrates.size() != M.size())
                    {
                        throw std::invalid_argument(
                            "migration matrix size mismatch");
                    }
                if (source != NULLDEME)
                    {
                        std::copy(begin(migrates), end(migrates),
                                  M.begin() + npops * source);
                    }
                else
                    {
                        std::copy(begin(migrates), end(migrates), M.begin());
                    }
            }

            bool
            empty() const
            {
                return M.empty();
            }
        
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
