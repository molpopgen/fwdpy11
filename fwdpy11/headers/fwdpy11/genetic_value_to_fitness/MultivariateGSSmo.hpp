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
#ifndef FWDPY11_MULTIVARIATE_GSSMO_HPP
#define FWDPY11_MULTIVARIATE_GSSMO_HPP

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "GeneticValueIsTrait.hpp"
#include "PleiotropicOptima.hpp"

namespace fwdpy11
{
    struct MultivariateGSSmo : public GeneticValueIsTrait
    {
        std::vector<PleiotropicOptima> optima;
        std::size_t current_timepoint;

        MultivariateGSSmo(const std::vector<PleiotropicOptima> &po)
            : GeneticValueIsTrait{po.empty() ? 0 : po[0].optima.size()}, optima(po),
              current_timepoint(0)
        {
            if (po.empty())
                {
                    throw std::invalid_argument("empty list of PleiotropicOptima");
                }
            for (auto &o : optima)
                {
                    if (o.when == PleiotropicOptima::null)
                        {
                            throw std::invalid_argument(
                                "invalid when value for PleiotropicOptima");
                        }
                }
            for (auto &o : optima)
                {
                    if (o.optima.size() != total_dim)
                        {
                            throw std::invalid_argument(
                                "all lists of optima must be the same length");
                        }
                }
            if (!std::is_sorted(begin(optima), end(optima)))
                {
                    throw std::invalid_argument("optima must be sorted by `when` field");
                }
        }

        double
        operator()(const DiploidGeneticValueToFitnessData data) const override
        {
            if (data.gvalues.get().size() != total_dim)
                {
                    throw std::runtime_error("dimension mismatch");
                }
            double sqdiff = 0.0;
            for (std::size_t i = 0; i < data.gvalues.get().size(); ++i)
                {
                    sqdiff += gsl_pow_2(data.gvalues.get()[i]
                                        - optima[current_timepoint].optima[i]);
                }
            return std::exp(-sqdiff / (2.0 * optima[current_timepoint].VW));
        }

        std::shared_ptr<GeneticValueToFitnessMap>
        clone() const override
        {
            return std::make_shared<MultivariateGSSmo>(optima);
        }

        template <typename poptype>
        inline void
        update_details(const poptype &pop)
        {
            if (current_timepoint < optima.size() - 1
                && pop.generation >= optima[current_timepoint].when)
                {
                    ++current_timepoint;
                }
        }

        void
        update(const DiploidPopulation &pop) override
        {
            update_details(pop);
        }
    };
} // namespace fwdpy11

#endif

