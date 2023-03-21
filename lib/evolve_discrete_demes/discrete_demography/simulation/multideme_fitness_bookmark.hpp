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

#include <system_error>
#include <vector>
#include <utility>
#include <numeric>
#include <limits>
#include <cstdint>
#include "core/demes/forward_graph.hpp"

namespace fwdpy11_core
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
            update(const fwdpy11_core::ForwardDemesGraphDataIterator<double> deme_sizes,
                   const std::vector<METADATATYPE>& individual_metadata)
            {
                std::vector<std::uint32_t> deme_sizes_uint;
                for (auto i = std::begin(deme_sizes); i != std::end(deme_sizes); ++i)
                    {
                        deme_sizes_uint.push_back(static_cast<std::uint32_t>(*i));
                    }
                auto ndemes = deme_sizes_uint.size();
                starts.resize(ndemes);
                stops.resize(ndemes);
                offsets.resize(ndemes);
                std::uint32_t ttl_N = std::accumulate(std::begin(deme_sizes_uint),
                                                      std::end(deme_sizes_uint), 0);
                individual_fitness.resize(ttl_N, -1.0);
                individuals.resize(individual_fitness.size(),
                                   std::numeric_limits<std::uint32_t>::max());
                // These -1 are useful b/c if our bookkeeping is bad,
                // we'll get errors from GSL.
                std::fill(begin(individual_fitness), end(individual_fitness), -1);
                std::partial_sum(std::begin(deme_sizes_uint), std::end(deme_sizes_uint),
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
