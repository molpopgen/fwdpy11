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
#ifndef FWDPY11_SLOCUSPOP_GENETIC_VALUE_WITH_MAPPING_HPP__
#define FWDPY11_SLOCUSPOP_GENETIC_VALUE_WITH_MAPPING_HPP__

#include <memory>
#include "SlocusPopGeneticValue.hpp"

namespace fwdpy11
{
    struct SlocusPopGeneticValueWithMapping : public SlocusPopGeneticValue
    ///API class.
    {
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::unique_ptr<GeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::unique_ptr<SlocusPopGeneticValueNoise> noise_fxn;

        SlocusPopGeneticValueWithMapping(
            std::unique_ptr<GeneticValueToFitnessMap> gv2w_)
            : gv2w{ std::move(gv2w_) }, noise_fxn{ new SlocusPopNoNoise() }
        {
        }

        SlocusPopGeneticValueWithMapping(
            std::unique_ptr<GeneticValueToFitnessMap> gv2w_,
            std::unique_ptr<SlocusPopGeneticValueNoise> noise_)
            : gv2w{ std::move(gv2w_) }, noise_fxn{ std::move(noise_) }
        {
        }

        inline virtual double
        genetic_value_to_fitness(const double g, const double e) const
        {
            return gv2w->operator()(g, e);
        }

        inline virtual double
        noise(const GSLrng_t& rng, const dip_metadata& offspring_metadata,
              const std::size_t parent1, const std::size_t parent2,
              const SlocusPop& pop) const
        {
            return noise_fxn->operator()(rng, offspring_metadata, parent1,
                                         parent2, pop);
        }
    };
} // namespace fwdpy11

#endif
