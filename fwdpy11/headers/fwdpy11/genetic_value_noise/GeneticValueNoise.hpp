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
#include <pybind11/pybind11.h>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct GeneticValueNoise
    ///ABC for random effects on trait values
    {
        virtual ~GeneticValueNoise() = default;
        virtual double
        operator()(const GSLrng_t& /* rng */,
                   const DiploidMetadata& /*offspring_metadata*/,
                   const std::size_t /*parent1*/,
                   const std::size_t /*parent2*/,
                   const DiploidPopulation& /*pop*/) const = 0;
        virtual void update(const DiploidPopulation& /*pop*/) = 0;
        virtual std::unique_ptr<GeneticValueNoise> clone() const = 0;
        virtual pybind11::object pickle() const = 0;
    };
} // namespace fwdpy11

#endif
