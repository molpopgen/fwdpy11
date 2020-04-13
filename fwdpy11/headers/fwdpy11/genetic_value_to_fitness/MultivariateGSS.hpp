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
#ifndef FWDPY11_MULTIVARIATE_GSS_HPP
#define FWDPY11_MULTIVARIATE_GSS_HPP

#include <gsl/gsl_math.h>
#include <cmath>
#include <vector>
#include "GeneticValueIsTrait.hpp"
#include "PleiotropicOptima.hpp"
#include "../genetic_values/default_update.hpp"

namespace fwdpy11
{
    struct MultivariateGSS : public GeneticValueIsTrait
    {
        std::vector<double> optima;
        double VS;

        MultivariateGSS(std::vector<double> input_optima, double VS_)
            : GeneticValueIsTrait{input_optima.size()}, optima(std::move(input_optima)),
              VS(VS_)
        {
        }

        explicit MultivariateGSS(const PleiotropicOptima& po)
            : GeneticValueIsTrait{po.optima.size()}, optima{po.optima}, VS{po.VW}
        {
        }

        virtual double
        operator()(const DiploidMetadata& /*metadata*/,
                   const std::vector<double>& values) const override
        {
            if (values.size() != optima.size())
                {
                    throw std::runtime_error("dimension mismatch");
                }
            double sqdiff = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i)
                {
                    sqdiff += gsl_pow_2(values[i] - optima[i]);
                }
            return std::exp(-sqdiff / (2.0 * VS));
        }

        std::unique_ptr<GeneticValueToFitnessMap>
        clone() const override
        {
            return std::unique_ptr<MultivariateGSS>(new MultivariateGSS(*this));
        }

        pybind11::object
        pickle() const override
        {
            pybind11::list l;
            for (auto x : optima)
                {
                    l.append(x);
                }
            return pybind11::make_tuple(l, VS);
        }

        DEFAULT_DIPLOID_POP_UPDATE()
    };
} // namespace fwdpy11

#endif

