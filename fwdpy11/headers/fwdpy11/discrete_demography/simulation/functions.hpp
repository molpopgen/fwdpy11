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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_FUNCTIONS_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_FUNCTIONS_HPP

#include <vector>
#include "../../rng.hpp"
#include "../current_event_state.hpp"
#include "../exceptions.hpp"
#include "core/demes/forward_graph.hpp"
#include "multideme_fitness_lookups.hpp"
#include "migration_lookup.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        namespace detail
        {
            inline void
            no_valid_parents(std::size_t i, std::uint32_t generation, std::uint32_t N)
            {
                std::ostringstream o;
                o << "deme " << i << " at time " << generation << " has size " << N
                  << " and no valid parents";
                throw DemographyError(o.str());
            }


        } // namespace detail

        template <typename METADATATYPE>
        inline void
        get_current_deme_sizes(const std::vector<METADATATYPE>& metadata,
                               current_deme_sizes_vector& deme_sizes)
        {
            auto& ref = deme_sizes.get();
            std::fill(begin(ref), end(ref), 0);
            for (auto&& i : metadata)
                {
                    ref[i.deme]++;
                }
        }

        inline void
        validate_parental_state(
            std::uint32_t generation,
            const multideme_fitness_lookups<std::uint32_t>& fitnesses,
            const fwdpy11_core::ForwardDemesGraph& demography)
        {
            if (demography.iterating_model())
                {
                    auto offspring_deme_sizes = demography.offspring_deme_sizes();
                    auto parental_deme_sizes = demography.parental_deme_sizes();

                    std::size_t offspring_deme = 0;
                    for (auto size = std::begin(offspring_deme_sizes);
                         size != std::end(offspring_deme_sizes);
                         ++size, ++offspring_deme)
                        {
                            if (*size > 0.0)
                                {
                                    // offspring deme exists
                                    auto ancestry_proportions
                                        = demography.offspring_ancestry_proportions(
                                            offspring_deme);
                                    std::uint32_t parental_deme = 0;
                                    for (auto prop = std::begin(ancestry_proportions);
                                         prop != std::end(ancestry_proportions);
                                         ++prop, ++parental_deme)
                                        {
                                            if (*prop > 0.0)
                                                {
                                                    auto parental_deme_size = *(
                                                        std::begin(parental_deme_sizes)
                                                        + parental_deme);
                                                    if (parental_deme_size == 0.0)
                                                        {
                                                            std::ostringstream o;
                                                            o << "deme "
                                                              << offspring_deme
                                                              << " has ancestry from "
                                                                 "deme "
                                                              << parental_deme
                                                              << " at time "
                                                              << generation
                                                              << " but the parental "
                                                                 "deme size "
                                                                 "is 0";
                                                            throw DemographyError(
                                                                o.str());
                                                        }
                                                    if (fitnesses.lookups[parental_deme]
                                                        == nullptr)
                                                        {
                                                            std::ostringstream o;
                                                            o << "fitness lookup table "
                                                                 "for "
                                                                 "parental deme "
                                                              << parental_deme
                                                              << " is empty at time "
                                                              << generation;
                                                            throw DemographyError(
                                                                o.str());
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }

    } // namespace discrete_demography
} // namespace fwdpy11

#endif
