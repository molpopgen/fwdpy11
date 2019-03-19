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
#ifndef FWDPY11_POP_GENETIC_VALUE_HPP__
#define FWDPY11_POP_GENETIC_VALUE_HPP__

#include <cstdint>
#include <vector>
#include <pybind11/pybind11.h>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include "GeneticValueToFitness.hpp"
#include "noise.hpp"

namespace fwdpy11
{
    struct DiploidPopulationGeneticValue
    /// API class
    /// For a single-locus simulation, we need the following concepts:
    /// 1. Calculate the genetic value of a diploid, g.  This is calculate_gvalue()
    /// 2. Calculate any random effects, or "noise", e.  This is noise(...)
    /// 3. Calculate the final fitness, w = f(g,e). This is genetic_value_to_fitness()
    /// 4. If a derived class has any of its own data, provide a means for updating. This is update()
    ///
    /// In a simulation, we want the following operations:
    /// pop.diploid_metadata[i].g = gv(i,pop);
    /// pop.diploid_metadata[i].e = gv.noise(rng, pop.diploid_medatadata[i],p1,p2,pop);
    /// pop.diploid_metadata[i].w = gv.genetic_value_to_fitness(pop.diploid_metadata[i]);
    /// These operations are handled by operator()
    /// At the end of a generation, update may be called.
    ///
    /// Things to note:
    /// Any deme/geography-specific details must be handled by the derived class.
    {
        std::size_t total_dim;
        mutable std::vector<double> gvalues;

        explicit DiploidPopulationGeneticValue(std::size_t dimensonality)
            : total_dim(dimensonality), gvalues(total_dim, 0.0)
        {
        }

        virtual ~DiploidPopulationGeneticValue() = default;

        // Callable from Python
        virtual double calculate_gvalue(const std::size_t /*diploid_index*/,
                                        const DiploidPopulation& /*pop*/) const = 0;
        // To be called from w/in a simulation
        virtual void
        operator()(const GSLrng_t& rng, std::size_t diploid_index,
                   const DiploidPopulation& pop, DiploidMetadata& metadata) const
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
                             const DiploidPopulation& /*pop*/) const = 0;
        virtual void update(const DiploidPopulation& /*pop*/) = 0;
        virtual pybind11::object pickle() const = 0;

        virtual pybind11::tuple shape() const = 0;

        std::vector<double>
        genetic_values() const
        {
            return gvalues;
        }
    };
} //namespace fwdpy11

#endif
