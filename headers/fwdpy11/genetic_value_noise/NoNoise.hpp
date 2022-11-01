//
// Copyright (C) 2029 Kevin Thornton <krthornt@uci.edu>
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
//
#ifndef FWDPY11_NO_NOISE_HPP
#define FWDPY11_NO_NOISE_HPP

#include "GeneticValueNoise.hpp"

namespace fwdpy11
{
    struct NoNoise : public GeneticValueNoise
    {
        double
        operator()(const DiploidGeneticValueNoiseData /*data*/) const override
        {
            return 0.;
        }

        void
        update(const DiploidPopulation& /*pop*/) override
        {
        }

        std::shared_ptr<GeneticValueNoise>
        clone() const override
        {
            return std::make_shared<NoNoise>();
        }
    };
} // namespace fwdpy11

#endif
