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
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpy11/rng.hpp>
#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "current_event_state.hpp"
#include "simulation/deme_properties.hpp"
#include "simulation/functions.hpp"
#include "simulation/multideme_fitness_lookups.hpp"
#include "simulation/migration_lookup.hpp"
#include "simulation/build_migration_lookup.hpp"
#include "simulation/multideme_fitness_bookmark.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename T>
        std::int32_t
        update_maxdeme_from_demography(std::int32_t m, const T& t)
        {
            for (auto&& i : t)
                {
                    m = std::max(m, i.deme);
                }
            return m;
        }

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
            MigrationMatrix M;
            std::int32_t maxdemes;
            deme_properties current_deme_parameters;
            multideme_fitness_bookmark fitness_bookmark;

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
                  maxdemes{-1}, current_deme_parameters{}, fitness_bookmark{}
            {
            }

            DiscreteDemographyState(const DiscreteDemographyState&) = default;
            DiscreteDemographyState& operator=(const DiscreteDemographyState&) = default;
            DiscreteDemographyState(DiscreteDemographyState&&) = default;

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

            // NOTE: this sets maxdemes
            void
            initialize(const fwdpy11::DiploidPopulation& pop)
            {
                std::int32_t maxdeme_from_metadata = -1;
                std::int32_t maxdeme_from_demography = -1;
                for (auto& md : pop.diploid_metadata)
                    {
                        if (md.deme < 0)
                            {
                                throw std::invalid_argument(
                                    "input deme labels must be non-negative");
                            }
                        maxdeme_from_metadata = std::max(maxdeme_from_metadata, md.deme);
                    }
                for (auto& m : mass_migrations.events)
                    {
                        maxdeme_from_demography
                            = std::max(maxdeme_from_demography, m.source);
                        maxdeme_from_demography
                            = std::max(maxdeme_from_demography, m.destination);
                    }
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, set_growth_rates.events);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, set_deme_sizes.events);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, set_selfing_rates.events);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, set_migration_rates.events);
                auto temp = std::max(maxdeme_from_metadata, maxdeme_from_demography) + 1;
                if (M.empty())
                    {
                        // no migration, so done
                        maxdemes = temp;
                    }
                else
                    {
                        if (static_cast<std::size_t>(temp) > M.npops)
                            {
                                throw std::invalid_argument(
                                    "MigrationMatrix contains too few demes");
                            }
                        if (static_cast<std::size_t>(temp) < M.npops)
                            {
                                throw std::invalid_argument(
                                    "MigrationMatrix contains too many demes");
                            }
                        maxdemes = std::max(temp, static_cast<std::int32_t>(M.npops));
                    }

                if (current_deme_parameters.current_deme_sizes.get().empty())
                    {
                        current_deme_parameters = deme_properties{
                            static_cast<std::uint32_t>(maxdemes), pop.diploid_metadata};
                    }
            }

            void
            early(const GSLrng_t& /*rng*/, std::uint32_t /*current_simulation_time*/,
                  std::vector<DiploidMetadata>& /*individual_metadata*/)
            {
            }

            void
            late(const GSLrng_t& rng, std::uint32_t current_simulation_time,
                 migration_lookup& miglookup,
                 std::vector<DiploidMetadata>& individual_metadata)
            {
                mass_migration(rng, current_simulation_time, mass_migrations,
                               current_deme_parameters.growth_rates,
                               current_deme_parameters.growth_rate_onset_times,
                               current_deme_parameters.growth_initial_sizes,
                               individual_metadata);
                get_current_deme_sizes(individual_metadata,
                                       current_deme_parameters.current_deme_sizes);
                fitness_bookmark.update(current_deme_parameters.current_deme_sizes,
                                        individual_metadata);
                std::copy(begin(current_deme_parameters.current_deme_sizes.get()),
                          end(current_deme_parameters.current_deme_sizes.get()),
                          begin(current_deme_parameters.next_deme_sizes.get()));
                detail::set_next_deme_sizes(current_simulation_time, set_deme_sizes,
                                            current_deme_parameters);
                detail::update_growth_rates(current_simulation_time, set_growth_rates,
                                            current_deme_parameters);
                detail::update_selfing_rates(current_simulation_time, set_selfing_rates,
                                             current_deme_parameters);
                detail::update_migration_matrix(current_simulation_time,
                                                set_migration_rates, M);
                auto next_global_N_ = detail::apply_growth_rates_get_next_global_N(
                    current_simulation_time, current_deme_parameters);
                set_next_global_N(next_global_N_);
                build_migration_lookup(M, current_deme_parameters.current_deme_sizes,
                                       current_deme_parameters.next_deme_sizes,
                                       miglookup);
            }
        };
    }
}

#endif
