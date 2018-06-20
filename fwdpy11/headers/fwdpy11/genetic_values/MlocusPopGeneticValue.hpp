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
#ifndef FWDPY11_MLOCUSPOP_GENETIC_VALUE_HPP__
#define FWDPY11_MLOCUSPOP_GENETIC_VALUE_HPP__

#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

namespace fwdpy11
{
    struct MlocusPopGeneticValue
    /// API class
    /// Concepts we have to deal with:
    /// Noise and aggregation are a bit trickier here?
    {
        virtual double operator()(const std::size_t /*diploid*/,
                                  const fwdpy11::MlocusPop& /*pop*/) const = 0;
        virtual double genetic_value_to_fitness(
            const DiploidMetadata& /*metadata*/) const = 0;
        virtual double noise(const GSLrng_t& /*rng*/,
                             const DiploidMetadata& /*offspring_metadata*/,
                             const std::size_t /*parent1*/,
                             const std::size_t /*parent2*/,
                             const MlocusPop& /*pop*/) const = 0;
        virtual void update(const fwdpy11::MlocusPop& /*pop*/) = 0;
    };
} // namespace fwdpy11

#endif
