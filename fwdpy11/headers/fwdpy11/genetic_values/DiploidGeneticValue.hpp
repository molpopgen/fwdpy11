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
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsFitness.hpp>
#include <fwdpy11/genetic_value_noise/NoNoise.hpp>
#include <fwdpy11/genetic_value_data/genetic_value_data.hpp>
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    struct DiploidGeneticValue
    /// API class
    ///
    /// Things to note:
    /// Any deme/geography-specific details must be handled by the derived class.
    {
        std::size_t total_dim;
        std::vector<double> gvalues;
        /// Classes deriving from this must call gv2w->update
        /// from their own update functions.
        std::shared_ptr<GeneticValueToFitnessMap> gv2w;
        /// This must be updated, too:
        std::shared_ptr<GeneticValueNoise> noise_fxn;

        explicit DiploidGeneticValue(std::size_t ndim)
            : total_dim(ndim), gvalues(total_dim, 0.),
              gv2w{new GeneticValueIsFitness{total_dim}}, noise_fxn{new NoNoise()}
        {
        }

        DiploidGeneticValue(std::size_t dimensonality,
                            const GeneticValueToFitnessMap& gv2w_)
            : total_dim(dimensonality),
              gvalues(total_dim, 0.0), gv2w{gv2w_.clone()}, noise_fxn{new NoNoise}
        {
        }

        DiploidGeneticValue(std::size_t dimensonality,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : total_dim(dimensonality),
              gvalues(total_dim, 0.0), gv2w{gv2w_.clone()}, noise_fxn{noise_.clone()}
        {
        }

        // The type is move-only
        virtual ~DiploidGeneticValue() = default;
        DiploidGeneticValue(const DiploidGeneticValue&) = delete;
        DiploidGeneticValue(DiploidGeneticValue&&) = default;
        DiploidGeneticValue& operator=(const DiploidGeneticValue&) = delete;

        virtual double calculate_gvalue(const DiploidGeneticValueData data) = 0;

        virtual void update(const DiploidPopulation& pop) = 0;

        // To be called from w/in a simulation
        inline void
        operator()(DiploidGeneticValueData data)
        {
            data.offspring_metadata.get().g = calculate_gvalue(data);
            data.offspring_metadata.get().e = noise(DiploidGeneticValueNoiseData(data));
            data.offspring_metadata.get().w = genetic_value_to_fitness(
                DiploidGeneticValueToFitnessData(data, gvalues));
        }

        virtual double
        genetic_value_to_fitness(const DiploidGeneticValueToFitnessData data)
        {
            return gv2w->operator()(data);
        }

        virtual double
        noise(const DiploidGeneticValueNoiseData data) const
        {
            return noise_fxn->operator()(data);
        }

        virtual pybind11::tuple
        shape() const
        {
            if (total_dim > 1 && total_dim != gvalues.size())
                {
                    throw std::runtime_error("dimensionality mismatch");
                }
            return pybind11::make_tuple(total_dim);
        }
    };
} //namespace fwdpy11

#endif
