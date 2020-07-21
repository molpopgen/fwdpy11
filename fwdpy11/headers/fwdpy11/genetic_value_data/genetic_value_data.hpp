//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWDPY11_GENETIC_DATA_HPP__
#define FWDPY11_GENETIC_DATA_HPP__

#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace fwdpy11
{
    struct DiploidGeneticValueData
    {
        std::reference_wrapper<const fwdpy11::GSLrng_t> rng;
        std::reference_wrapper<const fwdpy11::DiploidPopulation> pop;
        std::reference_wrapper<const fwdpy11::DiploidMetadata> parent1_metadata,
            parent2_metadata;
        std::size_t metadata_index;
        std::reference_wrapper<fwdpy11::DiploidMetadata> offspring_metadata;

        DiploidGeneticValueData(const fwdpy11::GSLrng_t& rng_,
                                const fwdpy11::DiploidPopulation& pop_,
                                const fwdpy11::DiploidMetadata& parent1_metadata_,
                                const fwdpy11::DiploidMetadata& parent2_metadata_,
                                std::size_t metadata_index_,
                                fwdpy11::DiploidMetadata& offspring_metadata_)
            : rng(rng_), pop(pop_), parent1_metadata(parent1_metadata_),
              parent2_metadata(parent2_metadata_), metadata_index(metadata_index_),
              offspring_metadata(offspring_metadata_)
        {
        }
    };

    // NOTE: the next two structs may be collapsible into one?

    struct DiploidGeneticValueNoiseData
    {
        std::reference_wrapper<const fwdpy11::GSLrng_t> rng;
        std::reference_wrapper<const fwdpy11::DiploidMetadata> parent1_metadata,
            parent2_metadata, offspring_metadata;
        std::size_t metadata_index;

        explicit DiploidGeneticValueNoiseData(const DiploidGeneticValueData& data)
            : rng(data.rng), parent1_metadata(data.parent1_metadata),
              parent2_metadata(data.parent2_metadata),
              offspring_metadata(data.offspring_metadata),
              metadata_index(data.metadata_index)
        {
        }
    };

    struct DiploidGeneticValueToFitnessData
    {
        std::reference_wrapper<const fwdpy11::GSLrng_t> rng;
        std::reference_wrapper<const fwdpy11::DiploidMetadata> parent1_metadata,
            parent2_metadata, offspring_metadata;
        std::reference_wrapper<const std::vector<double>> gvalues;
        std::size_t metadata_index;

        explicit DiploidGeneticValueToFitnessData(const DiploidGeneticValueData& data,
                                                  const std::vector<double>& gvalues_)
            : rng(data.rng), parent1_metadata(data.parent1_metadata),
              parent2_metadata(data.parent2_metadata),
              offspring_metadata(data.offspring_metadata),
              gvalues(gvalues_),
              metadata_index(data.metadata_index)
        {
        }
    };
}

#endif
