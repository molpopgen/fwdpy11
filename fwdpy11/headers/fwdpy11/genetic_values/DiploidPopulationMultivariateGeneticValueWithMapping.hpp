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
#ifndef FWDPY11_POP_MULTIVARIATE_GENETIC_VALUE_WITH_MAPPING_HPP__
#define FWDPY11_POP_MULTIVARIATE_GENETIC_VALUE_WITH_MAPPING_HPP__

#include <memory>
#include "DiploidPopulationGeneticValue.hpp"
#include "MultivariateGeneticValueToFitnessMap.hpp"
#include "noise.hpp"

namespace fwdpy11
{
    struct DiploidPopulationMultivariateGeneticValueWithMapping
        : public DiploidPopulationGeneticValue
    {
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<MultivariateGeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<GeneticValueNoise> noise_fxn;

        DiploidPopulationMultivariateGeneticValueWithMapping(
            std::size_t ndim,
            const MultivariateGeneticValueToFitnessMap& gv2w_)
            : DiploidPopulationGeneticValue(ndim), gv2w{ gv2w_.clone() }, noise_fxn{
                  new NoNoise()
              }
        {
        }

        DiploidPopulationMultivariateGeneticValueWithMapping(
            std::size_t ndim,
            const MultivariateGeneticValueToFitnessMap& gv2w_,
            const GeneticValueNoise& noise_)
            : DiploidPopulationGeneticValue(ndim), gv2w{ gv2w_.clone() }, noise_fxn{
                  noise_.clone()
              }
        {
        }

        DiploidPopulationMultivariateGeneticValueWithMapping(
            const DiploidPopulationMultivariateGeneticValueWithMapping& other)
            : DiploidPopulationGeneticValue(other.total_dim), gv2w{ other.gv2w->clone() },
              noise_fxn{ other.noise_fxn->clone() }
        {
        }

        virtual ~DiploidPopulationMultivariateGeneticValueWithMapping() = default;

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

        inline void
        update(const fwdpy11::DiploidPopulation& pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }

        pybind11::tuple
        shape() const
        {
            return pybind11::make_tuple(total_dim);
        }
    };
} // namespace fwdpy11

#endif

