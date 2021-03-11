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
#ifndef FWDPY11_GENETIC_VALUE_TO_FITNESS_MAP_HPP__
#define FWDPY11_GENETIC_VALUE_TO_FITNESS_MAP_HPP__

#include <memory>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>
#include <fwdpp/util/named_type.hpp>
#include <fwdpy11/genetic_value_data/genetic_value_data.hpp>

namespace fwdpy11
{
    struct genetic_value_maps_to_fitness
    {
    };

    using maps_to_fitness
        = fwdpp::strong_types::named_type<bool, genetic_value_maps_to_fitness>;

    struct GeneticValueToFitnessMap
    {
        std::size_t total_dim;
        bool isfitness;
        explicit GeneticValueToFitnessMap(std::size_t ndim, const maps_to_fitness& m)
            : total_dim{ndim}, isfitness{m.get()}
        {
        }

        virtual ~GeneticValueToFitnessMap() = default;
        GeneticValueToFitnessMap(const GeneticValueToFitnessMap&) = delete;
        GeneticValueToFitnessMap(GeneticValueToFitnessMap&&) = default;
        GeneticValueToFitnessMap& operator=(const GeneticValueToFitnessMap&) = delete;
        GeneticValueToFitnessMap& operator=(GeneticValueToFitnessMap&&) = default;

        virtual double
        operator()(const DiploidGeneticValueToFitnessData /*data*/) const = 0;
        virtual void update(const DiploidPopulation& /*pop*/) = 0;
        virtual std::shared_ptr<GeneticValueToFitnessMap> clone() const = 0;
    };
} //namespace fwdpy11

#endif
