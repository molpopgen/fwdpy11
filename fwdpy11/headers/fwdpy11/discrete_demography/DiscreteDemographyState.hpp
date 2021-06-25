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
#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
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
            template <typename T> struct events_with_range
            {
                std::vector<T> events;
                std::pair<std::size_t, std::size_t> event_range;
                template <typename Input>
                events_with_range(Input&& input)
                    : events(std::forward<Input>(input)), event_range{0, events.size()}
                {
                }
            };

            events_with_range<MassMigration> mass_migrations;
            events_with_range<SetExponentialGrowth> set_growth_rates;
            events_with_range<SetDemeSize> set_deme_sizes;
            events_with_range<SetSelfingRate> set_selfing_rates;
            MigrationMatrix M;
            events_with_range<SetMigrationRates> set_migration_rates;

            std::uint32_t next_global_N;

          public:
            const std::int32_t maxdemes;
            multideme_fitness_lookups<std::uint32_t> fitnesses;
            deme_properties current_deme_parameters;
            MigrationMatrix M;
            migration_lookup miglookup;

            // NOTE: demography.update_event_times() needs to have been
            // called first!
            template <typename METADATATYPE>
            DiscreteDemographyState(std::vector<METADATATYPE>& metadata,
                                    std::vector<MassMigration> mass_migrations,
                                    std::vector<SetExponentialGrowth> set_growth_rates,
                                    std::vector<SetDemeSize> size_changes,
                                    std::vector<SetSelfingRate> set_selfing_rates,
                                    MigrationMatrix M,
                                    std::vector<SetMigrationRates> set_migration_rates)
                : mass_migrations(std::move(mass_migrations)),
                  set_growth_rates{std::move(set_growth_rates)},
                  set_deme_sizes{std::move(set_deme_sizes)},
                  set_selfing_rates{std::move(set_selfing_rates)}, M{std::move(M)},
                  set_migration_rates{std::move(set_migration_rates)}, next_global_N(0),
                  maxdemes(get_max_number_of_demes()(metadata, demography)),
                  fitnesses(maxdemes), current_deme_parameters(maxdemes, metadata),
                  M(demography.migmatrix), miglookup(maxdemes, M.empty())
            {
            }

            DiscreteDemographyState(const DiscreteDemographyState& other)
                : next_global_N{other.next_global_N}, maxdemes{other.maxdemes},
                  fitnesses{other.fitnesses},
                  current_deme_parameters{other.current_deme_parameters}, M{other.M},
                  miglookup{maxdemes, M.empty()}
            {
                build_migration_lookup(M, current_deme_parameters.current_deme_sizes,
                                       miglookup);
            }

            // This constructor is only used when resetting
            // the state from an event like pickling a DiscreteDemography
            // instance.
            DiscreteDemographyState(std::int32_t maxdemes_,
                                    deme_properties current_deme_parameters_,
                                    MigrationMatrix M_)
                : next_global_N(0), maxdemes(maxdemes_), fitnesses(maxdemes),
                  current_deme_parameters(std::move(current_deme_parameters_)),
                  M(std::move(M_)), miglookup(maxdemes, M.empty())
            {
                next_global_N = std::accumulate(
                    begin(current_deme_parameters.next_deme_sizes.get()),
                    end(current_deme_parameters.next_deme_sizes.get()), 0u);
            }

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

            std::unique_ptr<demographic_model_state>
            clone() const
            {
                return std::make_unique<demographic_model_state>(*this);
            }
        };
    }
}

#endif
