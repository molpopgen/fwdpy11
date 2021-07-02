//
// Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_STATE_HPP
#define FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_STATE_HPP

#include <cstdint>
#include <vector>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "current_event_state.hpp"
#include "simulation/migration_lookup.hpp"
#include "simulation/multideme_fitness_lookups.hpp"
#include "simulation/deme_properties.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class DiscreteDemographyState
        /// Added in 0.16.0 to hold and manage
        /// the relevant data structures.
        {
          private:
            std::uint32_t next_global_N;

          public:
            current_event_state<MassMigration> mass_migrations;
            current_event_state<SetExponentialGrowth> set_growth_rates;
            current_event_state<SetDemeSize> set_deme_sizes;
            current_event_state<SetSelfingRate> set_selfing_rates;
            current_event_state<SetMigrationRates> set_migration_rates;
            // this is the input matrix.
            // The current state of the matrix at time "t"
            // is stored in "miglookup" below
            MigrationMatrix M;
            std::int32_t maxdemes;
            multideme_fitness_lookups<std::uint32_t> fitnesses;
            //deme_properties current_deme_parameters;
            migration_lookup_v2 miglookup;

            DiscreteDemographyState(std::vector<MassMigration> mass_migrations,
                                    std::vector<SetExponentialGrowth> set_growth_rates,
                                    std::vector<SetDemeSize> set_deme_sizes,
                                    std::vector<SetSelfingRate> set_selfing_rates,
                                    MigrationMatrix M,
                                    std::vector<SetMigrationRates> set_migration_rates)
                : next_global_N{0}, mass_migrations(std::move(mass_migrations)),
                  set_growth_rates{std::move(set_growth_rates)},
                  set_deme_sizes{std::move(set_deme_sizes)}, set_selfing_rates{std::move(
                                                                 set_selfing_rates)},
                  set_migration_rates{std::move(set_migration_rates)}, M{std::move(M)},
                  maxdemes{0}, fitnesses{0}, miglookup{M}
            // current_deme_parameters(maxdemes, metadata), miglookup{M}
            {
            }

            DiscreteDemographyState(const DiscreteDemographyState&) = default;
            DiscreteDemographyState& operator=(const DiscreteDemographyState&) = default;
            DiscreteDemographyState(DiscreteDemographyState&&) = default;

            // This constructor is only used when resetting
            // the state from an event like pickling a DiscreteDemography
            // instance.
            // DiscreteDemographyState(std::int32_t maxdemes_,
            //                         deme_properties current_deme_parameters_,
            //                         MigrationMatrix M_)
            //     : next_global_N(0), maxdemes(maxdemes_), fitnesses(maxdemes),
            //       current_deme_parameters(std::move(current_deme_parameters_)),
            //       M(std::move(M_)), miglookup(maxdemes, M.empty())
            // {
            //     next_global_N = std::accumulate(
            //         begin(current_deme_parameters.next_deme_sizes.get()),
            //         end(current_deme_parameters.next_deme_sizes.get()), 0u);
            // }

            void
            set_next_global_N(std::uint32_t N)
            {
                next_global_N = N;
            }

            std::uint32_t
            ttlN_next() const
            {
                return next_global_N;
            }

            bool
            will_go_globally_extinct() const
            {
                return ttlN_next() == 0;
            }

            void
            initialize(const DiploidPopulation&)
            {
                throw std::runtime_error("not implemented");
            }

            void
            early(const std::uint32_t /*current_simulation_time*/)
            {
            }

            void
            late(const std::uint32_t /*current_simulation_time*/)
            {
            }
        };
    }
}

#endif
