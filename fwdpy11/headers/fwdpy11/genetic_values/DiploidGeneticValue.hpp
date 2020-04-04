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
#ifndef FWDPY11_DIPLOID_GENETIC_VALUE_HPP__
#define FWDPY11_DIPLOID_GENETIC_VALUE_HPP__

#include <cstdint>
#include <vector>
#include <pybind11/pybind11.h>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsFitness.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>
#include <fwdpy11/genetic_value_noise/NoNoise.hpp>

namespace fwdpy11
{
    struct DiploidGeneticValue
    /// API class
    /// For a diploid simulation, we need the following concepts:
    /// 1. Calculate the genetic value of a diploid, g.  This is calculate_gvalue()
    /// 2. Calculate any random effects, or "noise", e.  This is noise(...)
    /// 3. Calculate the final fitness, w = f(g,e). This is genetic_value_to_fitness()
    /// 4. If a derived class has any of its own data, provide a means for updating. This is update()
    ///
    /// In a simulation, we want the following operations:
    /// pop.diploid_metadata[i].g = gv(i,metadata,pop);
    /// pop.diploid_metadata[i].e = gv.noise(rng, pop.diploid_medatadata[i],p1,p2,pop);
    /// pop.diploid_metadata[i].w = gv.genetic_value_to_fitness(pop.diploid_metadata[i], gv.gvalues);
    /// These operations are handled by operator()
    /// At the end of a generation, update may be called.
    ///
    /// Things to note:
    /// Any deme/geography-specific details must be handled by the derived class.
    {
        std::size_t total_dim;
        mutable std::vector<double> gvalues;
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<GeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<GeneticValueNoise> noise_fxn;

        explicit DiploidGeneticValue(std::size_t ndim)
            : total_dim(ndim), gvalues(total_dim, 0.),
              gv2w{ new GeneticValueIsFitness{ total_dim } }, noise_fxn{
                  new NoNoise()
              }
        {
        }

        DiploidGeneticValue(std::size_t dimensonality,
                            const GeneticValueToFitnessMap& gv2w_)
            : total_dim(dimensonality),
              gvalues(total_dim, 0.0), gv2w{ gv2w_.clone() }, noise_fxn{
                  new NoNoise
              }
        {
        }

        DiploidGeneticValue(std::size_t dimensonality,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : total_dim(dimensonality),
              gvalues(total_dim, 0.0), gv2w{ gv2w_.clone() }, noise_fxn{
                  noise_.clone()
              }
        {
        }

        virtual ~DiploidGeneticValue() = default;

        DiploidGeneticValue(const DiploidGeneticValue&) = delete;
        DiploidGeneticValue(DiploidGeneticValue&&) = default;
        DiploidGeneticValue& operator=(const DiploidGeneticValue&) = delete;

        // Callable from Python
        virtual double
        calculate_gvalue(const std::size_t /*diploid_index*/,
                         const DiploidMetadata & /*diploid_metadata*/,
                         const DiploidPopulation& /*pop*/) const = 0;
        virtual pybind11::object pickle() const = 0;

        virtual void
        update(const DiploidPopulation& pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }

        // To be called from w/in a simulation
        virtual void
        operator()(const GSLrng_t& rng, std::size_t diploid_index,
                   const DiploidPopulation& pop,
                   DiploidMetadata& metadata) const
        {
            metadata.g = calculate_gvalue(diploid_index, metadata, pop);
            metadata.e = noise(rng, metadata, metadata.parents[0],
                               metadata.parents[1], pop);
            metadata.w = genetic_value_to_fitness(metadata);
        }

        virtual double
        genetic_value_to_fitness(const DiploidMetadata& metadata) const
        {
            return gv2w->operator()(metadata, gvalues);
        }

        virtual double
        noise(const GSLrng_t& rng, const DiploidMetadata& offspring_metadata,
              const std::size_t parent1, const std::size_t parent2,
              const DiploidPopulation& pop) const
        {
            return noise_fxn->operator()(rng, offspring_metadata, parent1,
                                         parent2, pop);
        }

        virtual pybind11::tuple
        shape() const
        {
            if (total_dim>1 && total_dim != gvalues.size())
                {
                    throw std::runtime_error("dimensionality mismatch");
                }
            return pybind11::make_tuple(total_dim);
        }
    };
} //namespace fwdpy11

#endif
