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
#include "MultivariateGeneticValueToFitnessMap.hpp"

namespace fwdpy11
{
    struct MultivariateGSSmo : public MultivariateGeneticValueToFitnessMap
    {
        std::vector<std::uint32_t> timepoints;
        std::vector<double> optima;
        std::size_t current_timepoint, ndim, optima_offset;
        double VS;

        MultivariateGSSmo(std::vector<std::uint32_t> input_timepoints,
                          std::vector<double> input_optima, double VS_)
            : timepoints(std::move(input_timepoints)),
              optima(std::move(input_optima)), current_timepoint(1), ndim(0), optima_offset(0),
              VS(VS_)
        {
            if (timepoints.empty())
                {
                    throw std::invalid_argument("empty timepoints");
                }
            if (optima.empty())
                {
                    throw std::invalid_argument("empty optima");
                }
            if (timepoints.front() != 0)
                {
                    throw std::invalid_argument(
                        "first timepoint is not at zero");
                }
            if (optima.size() % timepoints.size() != 0.0)
                {
                    throw std::invalid_argument(
                        "incorrect number of optima or time points");
                }
            if (!std::isfinite(VS))
                {
                    throw std::invalid_argument("VS must be finite");
                }
            if (VS <= 0.0)
                {
                    throw std::invalid_argument("VS must be >= 0");
                }
            ndim = optima.size() / timepoints.size();
        }

        virtual double
        operator()(const DiploidMetadata & /*metadata*/,
                   const std::vector<double> &values) const
        {
            if (values.size() != ndim)
                {
                    throw std::runtime_error("dimension mismatch");
                }
            double sqdiff = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i)
                {
                    sqdiff += gsl_pow_2(
                        values[i] - optima[optima_offset + i]);
                }
            return std::exp(-sqdiff / (2.0 * VS));
        }

        std::unique_ptr<MultivariateGeneticValueToFitnessMap>
        clone() const
        {
            return std::unique_ptr<MultivariateGSSmo>(
                new MultivariateGSSmo(*this));
        }

        pybind11::object
        pickle() const
        {
            pybind11::list l;
            for (auto x : optima)
                {
                    l.append(x);
                }
            pybind11::list tp;
            for (auto x : timepoints)
                {
                    tp.append(x);
                }
            return pybind11::make_tuple(tp, l, VS);
        }

        template <typename poptype>
        inline void
        update_details(const poptype &pop)
        {
            if (current_timepoint < timepoints.size()
                && pop.generation >= timepoints[current_timepoint])
                {
                    optima_offset += ndim;
                    ++current_timepoint;
                }
        }

        inline void
        update(const DiploidPopulation &pop)
        {
            update_details(pop);
        }
    };
} // namespace fwdpy11

#endif

