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
#ifndef FWDPY11_POP_GENETIC_VALUE_WITH_MAPPING_HPP__
#define FWDPY11_POP_GENETIC_VALUE_WITH_MAPPING_HPP__

#include <memory>
#include "DiploidPopulationGeneticValue.hpp"

namespace fwdpy11
{
    struct DiploidPopulationGeneticValueWithMapping : public DiploidPopulationGeneticValue
    ///API class.
    {
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<GeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<GeneticValueNoise> noise_fxn;

        DiploidPopulationGeneticValueWithMapping(const GeneticValueToFitnessMap& gv2w_)
            : DiploidPopulationGeneticValue(1), gv2w{ gv2w_.clone() },
              noise_fxn{ new NoNoise() }
        {
        }

        DiploidPopulationGeneticValueWithMapping(const GeneticValueToFitnessMap& gv2w_,
                                         const GeneticValueNoise& noise_)
            : DiploidPopulationGeneticValue(1), gv2w{ gv2w_.clone() }, noise_fxn{
                  noise_.clone()
              }
        {
        }

        inline virtual double
        genetic_value_to_fitness(const DiploidMetadata& metadata) const
        {
            return gv2w->operator()(metadata);
        }

        inline virtual double
        noise(const GSLrng_t& rng, const DiploidMetadata& offspring_metadata,
              const std::size_t parent1, const std::size_t parent2,
              const DiploidPopulation& pop) const
        {
            return noise_fxn->operator()(rng, offspring_metadata, parent1,
                                         parent2, pop);
        }

        pybind11::tuple
        shape() const
        {
            if (total_dim != 1 || total_dim != gvalues.size())
                {
                    throw std::runtime_error("dimensionality mismatch");
                }
            return pybind11::make_tuple(total_dim);
        }
    };
} // namespace fwdpy11

#endif
