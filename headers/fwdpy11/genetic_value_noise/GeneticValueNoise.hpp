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
#ifndef FWDPY11_GENETIC_VALUES_NOISE_HPP__
#define FWDPY11_GENETIC_VALUES_NOISE_HPP__

#include <memory>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/genetic_value_data/genetic_value_data.hpp>

namespace fwdpy11
{
    struct GeneticValueNoise
    ///ABC for random effects on trait values
    {
        virtual ~GeneticValueNoise() = default;
        GeneticValueNoise() = default;
        GeneticValueNoise(const GeneticValueNoise &) = delete;
        GeneticValueNoise(GeneticValueNoise &&) = default;
        GeneticValueNoise& operator=(GeneticValueNoise &&) = default;
        GeneticValueNoise& operator=(const GeneticValueNoise &) = delete;
        virtual double
        operator()(const DiploidGeneticValueNoiseData /*data*/) const = 0;
        virtual void update(const DiploidPopulation& /*pop*/) = 0;
        virtual std::shared_ptr<GeneticValueNoise> clone() const = 0;
    };
} // namespace fwdpy11

#endif
