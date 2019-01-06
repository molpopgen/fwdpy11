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
#ifndef FWDPY11_SLOCUSPOP_MULTIVARIATE_GENETIC_VALUE_WITH_MAPPING_HPP__
#define FWDPY11_SLOCUSPOP_MULTIVARIATE_GENETIC_VALUE_WITH_MAPPING_HPP__

#include <memory>
#include "SlocusPopGeneticValue.hpp"
#include "MultivariateGeneticValueToFitnessMap.hpp"
#include "noise.hpp"

namespace fwdpy11
{
    struct SlocusPopMultivariateGeneticValueWithMapping
        : public SlocusPopGeneticValue
    {
        mutable std::vector<double> multivariate_effects;
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<MultivariateGeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<GeneticValueNoise> noise_fxn;

        SlocusPopMultivariateGeneticValueWithMapping(
            std::size_t ndim,
            const MultivariateGeneticValueToFitnessMap& gv2w_)
            : multivariate_effects(ndim, 0.0), gv2w{ gv2w_.clone() },
              noise_fxn{ new NoNoise() }
        {
        }

        SlocusPopMultivariateGeneticValueWithMapping(
            std::size_t ndim,
            const MultivariateGeneticValueToFitnessMap& gv2w_,
            const GeneticValueNoise& noise_)
            : multivariate_effects(ndim, 0.0), gv2w{ gv2w_.clone() },
              noise_fxn{ noise_.clone() }
        {
        }

        virtual double
        genetic_value_to_fitness(const DiploidMetadata& metadata) const
        {
            return gv2w->operator()(metadata, multivariate_effects);
        }

        virtual double
        noise(const GSLrng_t& rng, const DiploidMetadata& offspring_metadata,
              const std::size_t parent1, const std::size_t parent2,
              const SlocusPop& pop) const
        {
            return noise_fxn->operator()(rng, offspring_metadata, parent1,
                                         parent2, pop);
        }

        inline void
        update(const fwdpy11::SlocusPop& pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }
    };
} // namespace fwdpy11

#endif

