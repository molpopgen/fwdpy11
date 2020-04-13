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

#ifndef FWDPY11_GSS
#define FWDPY11_GSS

#include <cmath>
#include <vector>
#include "GeneticValueIsTrait.hpp"
#include "Optimum.hpp"

namespace fwdpy11
{
    struct GSS : public GeneticValueIsTrait
    {
        const double opt, VS;
        GSS(const double opt_, const double VS_)
            : GeneticValueIsTrait{1}, opt{opt_}, VS{VS_}
        {
            if (VS <= 0.0)
                {
                    throw std::invalid_argument("VS must be > 0.0");
                }
            if (!std::isfinite(VS) || !std::isfinite(opt))
                {
                    throw std::invalid_argument("Both VS and opt must be finite values");
                }
        }

        explicit GSS(const Optimum &o) : GeneticValueIsTrait{1}, opt{o.opt}, VS{o.VW}
        {
        }

        double
        operator()(const DiploidMetadata &metadata,
                   const std::vector<double> & /*genetic_values*/) const override
        {
            return std::exp(
                -(std::pow(metadata.g + metadata.e - opt, 2.0) / (2.0 * VS)));
        }

        void
        update(const DiploidPopulation & /*pop*/) override
        {
        }

        std::unique_ptr<GeneticValueToFitnessMap>
        clone() const override
        {
            return std::unique_ptr<GSS>(new GSS(opt, VS));
        }

        pybind11::object
        pickle() const override
        {
            return pybind11::make_tuple(opt, VS);
        }
    };
} // namespace fwdpy11

#endif
