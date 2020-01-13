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

#ifndef FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_MANAGER_HPP
#define FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_MANAGER_HPP

#include <memory>
#include <cstdint>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include "get_max_number_of_demes.hpp"
#include "multideme_fitness_lookups.hpp"
#include "migration_lookup.hpp"
#include "deme_properties.hpp"
#include "../DiscreteDemography.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class discrete_demography_manager
        /// Added in 0.6.0 to hold and manage
        /// the relevant data structures.
        /// A current limitation is that only 
        /// DiploidPopulation is supported.
        {
          private:
            std::unique_ptr<MigrationMatrix>
            init_migmatrix(
                const std::unique_ptr<const MigrationMatrix> &Minput)
            {
                if (Minput == nullptr)
                    {
                        return nullptr;
                    }
                return std::unique_ptr<MigrationMatrix>(
                    new MigrationMatrix(*Minput));
            }

          public:
            const std::int32_t maxdemes;
            multideme_fitness_lookups<std::uint32_t> fitnesses;
            deme_properties sizes_rates;
            std::unique_ptr<MigrationMatrix> M;
            migration_lookup miglookup;

            // NOTE: demography.update_event_times() needs to have been
            // called first!
            discrete_demography_manager(const fwdpy11::DiploidPopulation &pop,
                                        DiscreteDemography &demography)
                : maxdemes(get_max_number_of_demes()(pop.diploid_metadata,
                                                     demography)),
                  fitnesses(maxdemes),
                  sizes_rates(maxdemes, pop.diploid_metadata),
                  M(init_migmatrix(demography.migmatrix)),
                  miglookup(maxdemes, M == nullptr)
            {
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
