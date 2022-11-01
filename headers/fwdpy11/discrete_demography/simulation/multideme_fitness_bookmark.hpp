//
// Copyright (C) 2021 is part of fwdpy11.
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

#pragma once

#include <vector>
#include <utility>
#include <numeric>
#include <limits>
#include <cstdint>
#include "deme_property_types.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct multideme_fitness_bookmark
        {
            std::vector<std::uint32_t> starts, stops, offsets, individuals;
            std::vector<double> individual_fitness;

            multideme_fitness_bookmark()
                : starts{}, stops{}, offsets{}, individuals{}, individual_fitness{}
            {
            }

            multideme_fitness_bookmark(std::vector<std::uint32_t> starts,
                                       std::vector<std::uint32_t> stops,
                                       std::vector<std::uint32_t> offsets,
                                       std::vector<std::uint32_t> individuals,
                                       std::vector<double> individual_fitness)
                : starts{std::move(starts)}, stops{std::move(stops)}, offsets{std::move(
                                                                          offsets)},
                  individuals{std::move(individuals)}, individual_fitness{
                                                           std::move(individual_fitness)}
            {
            }

            template <typename METADATATYPE>
            void
            update(const current_deme_sizes_vector& deme_sizes,
                   const std::vector<METADATATYPE>& individual_metadata)
            {
                auto& deme_sizes_ref = deme_sizes.get();
                starts.resize(deme_sizes_ref.size());
                stops.resize(deme_sizes_ref.size());
                offsets.resize(deme_sizes_ref.size());
                individual_fitness.resize(
                    std::accumulate(begin(deme_sizes_ref), end(deme_sizes_ref), 0),
                    -1.0);
                individuals.resize(individual_fitness.size(),
                                   std::numeric_limits<std::uint32_t>::max());
                // These -1 are useful b/c if our bookkeeping is bad,
                // we'll get errors from GSL.
                std::fill(begin(individual_fitness), end(individual_fitness), -1);
                std::partial_sum(begin(deme_sizes_ref), end(deme_sizes_ref),
                                 begin(stops));
                std::copy(begin(stops), end(stops) - 1, begin(starts) + 1);
                std::fill(begin(offsets), end(offsets), 0);
                for (auto&& md : individual_metadata)
                    {
                        auto i = starts[md.deme] + offsets[md.deme];
                        individual_fitness[i] = md.w;
                        individuals[i] = md.label;
                        offsets[md.deme]++;
                    }
            }
        };
    }
}
