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
#ifndef FWDPY11_MULTIVARIATE_GENETIC_VALUE_TO_FITNESS_HPP__
#define FWDPY11_MULTIVARIATE_GENETIC_VALUE_TO_FITNESS_HPP__

#include <cmath>
#include <memory>
#include <functional>
#include <algorithm>
#include <limits>
#include <vector>
#include <queue>
#include <tuple>
#include <pybind11/pybind11.h>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace fwdpy11
{
    struct MultivariateGeneticValueToFitnessMap
    {
        virtual ~MultivariateGeneticValueToFitnessMap() = default;
        virtual double
        operator()(const DiploidMetadata& /*metadata*/,
                   const std::vector<double>& /*values*/) const = 0;
        virtual void update(const DiploidPopulation& /*pop*/) = 0;
        virtual std::unique_ptr<MultivariateGeneticValueToFitnessMap>
        clone() const = 0;
        virtual pybind11::object pickle() const = 0;
    };
} // namespace fwdpy11

#endif
