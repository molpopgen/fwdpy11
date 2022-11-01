//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWDPY11_GENETIC_VALUE_IS_TRAIT
#define FWDPY11_GENETIC_VALUE_IS_TRAIT

#include "GeneticValueToFitnessMap.hpp"

namespace fwdpy11
{
    struct GeneticValueIsTrait : public GeneticValueToFitnessMap
    /// Another ABC.  Effectively a type trait
    {
        explicit GeneticValueIsTrait(std::size_t ndim)
            : GeneticValueToFitnessMap(ndim, maps_to_fitness(false))
        {
        }
        virtual ~GeneticValueIsTrait() = default;
        GeneticValueIsTrait(const GeneticValueIsTrait&) = delete;
        GeneticValueIsTrait(GeneticValueIsTrait&&) = default;
        GeneticValueIsTrait& operator=(const GeneticValueIsTrait&) = delete;
        GeneticValueIsTrait& operator=(GeneticValueIsTrait&&) = default;
    };
} // namespace fwdpy11

#endif
