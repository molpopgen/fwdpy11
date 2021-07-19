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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_PICK_PARENTS_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_PICK_PARENTS_HPP

#include <cstdint>
#include <sstream>
#include "../../rng.hpp"
#include "../MigrationMatrix.hpp"
#include "migration_lookup.hpp"
#include "deme_property_types.hpp"
#include "multideme_fitness_lookups.hpp"
#include "mating_event_type.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct parent_data
        {
            std::size_t parent1, parent2;
            std::int32_t deme1, deme2;
            mating_event_type mating;

            parent_data(std::size_t p1, std::size_t p2, std::int32_t d1, std::int32_t d2,
                        mating_event_type m)
                : parent1(p1), parent2(p2), deme1(d1), deme2(d2), mating(m)
            {
            }
        };

        inline parent_data
        pick_parents(const GSLrng_t& rng, const std::int32_t offspring_deme,
                     const migration_lookup& miglookup,
                     const current_deme_sizes_vector& current_deme_sizes,
                     const selfing_rates_vector& selfing_rates,
                     const multideme_fitness_bookmark& fitness_bookmark,
                     const multideme_fitness_lookups<std::uint32_t>& wlookups)
        {
            if (miglookup.null_migmatrix) // Model has no migration
                {
                    auto p1 = wlookups.get_parent(rng, current_deme_sizes,
                                                  fitness_bookmark, offspring_deme);
                    if (selfing_rates.get()[offspring_deme] > 0.
                        && gsl_rng_uniform(rng.get())
                               <= selfing_rates.get()[offspring_deme])
                        {
                            return {p1, p1, offspring_deme, offspring_deme,
                                    mating_event_type::selfing};
                        }
                    auto p2 = wlookups.get_parent(rng, current_deme_sizes,
                                                  fitness_bookmark, offspring_deme);
                    return {p1, p2, offspring_deme, offspring_deme,
                            mating_event_type::outcrossing};
                }
            if (miglookup.lookups[offspring_deme] == nullptr)
                {
                    std::ostringstream o;
                    o << "parental deme lookup is NULL for deme " << offspring_deme;
                    throw DemographyError(o.str());
                }
            std::int32_t pdeme = static_cast<std::int32_t>(
                gsl_ran_discrete(rng.get(), miglookup.lookups[offspring_deme].get()));

            auto p1
                = wlookups.get_parent(rng, current_deme_sizes, fitness_bookmark, pdeme);
            if (selfing_rates.get()[pdeme] > 0.
                && gsl_rng_uniform(rng.get()) <= selfing_rates.get()[pdeme])
                {
                    return {p1, p1, pdeme, pdeme, mating_event_type::selfing};
                }
            auto p2
                = wlookups.get_parent(rng, current_deme_sizes, fitness_bookmark, pdeme);
            return {p1, p2, pdeme, pdeme, mating_event_type::outcrossing};
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
