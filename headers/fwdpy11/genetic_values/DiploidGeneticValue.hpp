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

namespace fwdpy11
{
    class DiploidGeneticValue
    /// API class
    ///
    /// Things to note:
    /// Any deme/geography-specific details must be handled by the derived class.
    {
      private:
        template <typename Base, typename Default, typename... Args>
        std::shared_ptr<Base>
        process_input(const Base* input, Args... args)
        {
            if (input == nullptr)
                {
                    return std::shared_ptr<Base>(new Default{args...});
                }
            return input->clone();
        }

      public:
        std::size_t total_dim;
        std::vector<double> gvalues;
        // Even though these are stored as shared_ptr,
        // this class is non-copyable because its state
        // may change over time via the various update
        // functions.
        std::shared_ptr<GeneticValueToFitnessMap> gv2w;
        std::shared_ptr<GeneticValueNoise> noise_fxn;

        DiploidGeneticValue(std::size_t ndim, const GeneticValueToFitnessMap* gv2w_,
                            const GeneticValueNoise* noise)
            : total_dim(ndim), gvalues(total_dim, 0.),
              gv2w{process_input<GeneticValueToFitnessMap, GeneticValueIsFitness,
                                 std::size_t>(gv2w_, ndim)},
              noise_fxn{process_input<GeneticValueNoise, NoNoise>(noise)}
        {
        }

        // The type is move-only
        virtual ~DiploidGeneticValue() = default;
        DiploidGeneticValue(const DiploidGeneticValue&) = delete;
        DiploidGeneticValue(DiploidGeneticValue&&) = default;
        DiploidGeneticValue& operator=(const DiploidGeneticValue&) = delete;
        DiploidGeneticValue& operator=(DiploidGeneticValue&&) = default;

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
    };
} //namespace fwdpy11

#endif
