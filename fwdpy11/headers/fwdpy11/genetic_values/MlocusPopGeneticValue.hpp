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
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    struct MlocusPopGeneticValue
    /// API class
    /// Concepts we have to deal with:
    /// Noise and aggregation are a bit trickier here?
    {
        std::size_t total_dim;
        mutable std::vector<double> gvalues;

        explicit MlocusPopGeneticValue(std::size_t dimensonality)
            : total_dim(dimensonality), gvalues(total_dim, 0.0)
        {
        }

        virtual ~MlocusPopGeneticValue() = default;

        // Callable from Python
        virtual double
        calculate_gvalue(const std::size_t /*diploid*/,
                         const fwdpy11::MlocusPop& /*pop*/) const = 0;

        // To be called from w/in a simulation
        virtual void
        operator()(const GSLrng_t& rng, std::size_t diploid_index,
                   const MlocusPop& pop, DiploidMetadata& metadata) const
        {
            metadata.g = calculate_gvalue(diploid_index, pop);
            metadata.e = noise(rng, metadata, metadata.parents[0],
                               metadata.parents[1], pop);
            metadata.w = genetic_value_to_fitness(metadata);
        }

        virtual double genetic_value_to_fitness(
            const DiploidMetadata& /*metadata*/) const = 0;
        virtual double noise(const GSLrng_t& /*rng*/,
                             const DiploidMetadata& /*offspring_metadata*/,
                             const std::size_t /*parent1*/,
                             const std::size_t /*parent2*/,
                             const MlocusPop& /*pop*/) const = 0;
        virtual void update(const fwdpy11::MlocusPop& /*pop*/) = 0;
        virtual pybind11::tuple shape() const = 0;

        std::vector<double>
        genetic_values() const
        {
            return gvalues;
        }
        virtual pybind11::object pickle() const = 0;
    };
} // namespace fwdpy11

#endif
