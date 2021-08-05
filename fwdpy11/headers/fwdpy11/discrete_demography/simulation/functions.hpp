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
#include "apply_mass_migrations.hpp"
#include "../MassMigration.hpp"
#include "../exceptions.hpp"
#include "multideme_fitness_lookups.hpp"
#include "migration_lookup.hpp"
#include "deme_properties.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename METADATATYPE>
        inline void
        mass_migration(const GSLrng_t& rng, std::uint32_t t,
                       current_event_state<MassMigration>& mass_migration_events,
                       growth_rates_vector& growth_rates,
                       growth_rates_onset_times_vector& growth_rate_onset_times,
                       growth_initial_size_vector& growth_initial_sizes,
                       std::vector<METADATATYPE>& metadata)
        {
            if (mass_migration_events.current() < mass_migration_events.last()
                && t < mass_migration_events.when())
                {
                    return;
                }
            while (mass_migration_events.current() < mass_migration_events.last()
                   && mass_migration_events.when() == t)
                {
                    apply_mass_migrations(rng, t, mass_migration_events, growth_rates,
                                          growth_rate_onset_times, growth_initial_sizes,
                                          metadata);
                }
        }

        namespace detail
        {
            inline void
            set_next_deme_sizes(const std::uint32_t t,
                                current_event_state<SetDemeSize>& size_change_events,
                                deme_properties& sizes_rates)
            // NOTE: this function may resent growth rates to zero.
            // see SetDemeSize for details.
            {
                if (size_change_events.current() < size_change_events.last()
                    && t < size_change_events.when())
                    {
                        return;
                    }
                next_deme_sizes_vector::value_type& next_deme_sizes
                    = sizes_rates.next_deme_sizes.get();
                growth_rates_vector::value_type& growth_rates
                    = sizes_rates.growth_rates.get();
                growth_rates_onset_times_vector::value_type& growth_rate_onset_times
                    = sizes_rates.growth_rate_onset_times.get();
                growth_initial_size_vector::value_type& growth_initial_sizes
                    = sizes_rates.growth_initial_sizes.get();
                for (; size_change_events.current() < size_change_events.last()
                       && size_change_events.when() == t;
                     ++size_change_events.current())
                    {
                        auto deme = size_change_events.event().deme;
                        next_deme_sizes[deme] = size_change_events.event().new_size;
                        if (size_change_events.event().resets_growth_rate == true)
                            {
                                growth_rates[deme] = NOGROWTH;
                            }
                        // Deme size has changed, so we reset
                        // the onset of growth and the initial N for this deme
                        // NOTE: this needs careful documentation!
                        growth_rate_onset_times[deme] = t;
                        growth_initial_sizes[deme] = next_deme_sizes[deme];
                    }
            }

            inline void
            update_growth_rates(
                const std::uint32_t t,
                current_event_state<SetExponentialGrowth>& growth_rate_changes,
                deme_properties& sizes_rates)
            {
                if (growth_rate_changes.current() < growth_rate_changes.last()
                    && t < growth_rate_changes.when())
                    {
                        return;
                    }
                auto& rates = sizes_rates.growth_rates.get();
                auto& onsets = sizes_rates.growth_rate_onset_times.get();
                auto& next_deme_sizes = sizes_rates.next_deme_sizes.get();
                auto& N0 = sizes_rates.growth_initial_sizes.get();
                for (; growth_rate_changes.current() < growth_rate_changes.last()
                       && growth_rate_changes.when() == t;
                     ++growth_rate_changes.current())
                    {
                        auto deme = growth_rate_changes.event().deme;
                        if (growth_rate_changes.event().G != NOGROWTH
                            && next_deme_sizes[deme] == 0)
                            {
                                throw DemographyError(
                                    "attempt to change growth rate in extinct "
                                    "deme");
                            }
                        rates[deme] = growth_rate_changes.event().G;
                        onsets[deme] = t;
                        N0[deme] = next_deme_sizes[deme];
                    }
            }

            inline void
            update_selfing_rates(
                const std::uint32_t t,
                current_event_state<SetSelfingRate>& selfing_rate_changes,
                deme_properties& sizes_rates)
            {
                if (selfing_rate_changes.current() < selfing_rate_changes.last()
                    && t < selfing_rate_changes.when())
                    {
                        return;
                    }
                auto& rates = sizes_rates.selfing_rates.get();
                auto& Nnext = sizes_rates.next_deme_sizes.get();
                for (; selfing_rate_changes.current() < selfing_rate_changes.last()
                       && selfing_rate_changes.when() == t;
                     ++selfing_rate_changes.current())
                    {
                        auto& event = selfing_rate_changes.event();
                        if (Nnext[event.deme] == 0)
                            {
                                throw DemographyError("attempt to set selfing "
                                                      "rate in extinct deme");
                            }
                        rates[event.deme] = event.S;
                    }
            }

            inline void
            update_migration_matrix(
                const std::uint32_t t,
                current_event_state<SetMigrationRates>& migration_rate_changes,
                MigrationMatrix& M)

            {
                if (migration_rate_changes.current() < migration_rate_changes.last()
                    && t < migration_rate_changes.when())
                    {
                        return;
                    }
                for (; migration_rate_changes.current() < migration_rate_changes.last()
                       && migration_rate_changes.when() == t;
                     ++migration_rate_changes.current())
                    {
                        auto& event = migration_rate_changes.event();
                        M.set_migration_rates(event.deme, event.migrates);
                    }
            }

            inline std::uint32_t
            apply_growth_rates_get_next_global_N(const std::uint32_t t,
                                                 deme_properties& sizes_rates)
            {
                auto& Nnext = sizes_rates.next_deme_sizes.get();
                auto& G = sizes_rates.growth_rates.get();
                auto& N0 = sizes_rates.growth_initial_sizes.get();
                auto& onset = sizes_rates.growth_rate_onset_times.get();
                std::uint32_t next_global_N = 0;
                for (std::size_t deme = 0; deme < Nnext.size(); ++deme)
                    {
                        if (G[deme] != NOGROWTH)
                            {
                                if (Nnext[deme] == 0)
                                    {
                                        throw DemographyError(
                                            "growth is happening in an "
                                            "extinct deme");
                                    }
                                double next_size = std::round(
                                    static_cast<double>(N0[deme])
                                    * std::pow(G[deme], static_cast<double>(
                                                            t - onset[deme] + 1)));
                                if (next_size <= 0.0)
                                    {
                                        next_size = 0.0;
                                        G[deme] = NOGROWTH;
                                        onset[deme] = t;
                                        // If the deme is recolonized,
                                        // this will be reset
                                        N0[deme] = 0;
                                    }
                                Nnext[deme] = next_size;
                            }
                        next_global_N += Nnext[deme];
                    }
                return next_global_N;
            }

            inline void
            no_valid_parents(std::size_t i, std::uint32_t generation, std::uint32_t N)
            {
                std::ostringstream o;
                o << "deme " << i << " at time " << generation << " has size " << N
                  << " and no valid parents";
                throw DemographyError(o.str());
            }

            inline void
            check_migration_in(std::size_t i, std::uint32_t generation,
                               const MigrationMatrix& M)
            {
                if (M.M[i * M.npops + i] > 0)
                    {
                        std::ostringstream o;
                        o << "deme " << i << " at time " << generation
                          << " has no valid parents "
                             "from "
                             "the same deme but "
                             "M[i,i] != "
                             "0";
                        throw DemographyError(o.str());
                    }
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
            const deme_properties& current_deme_parameters, const MigrationMatrix& M)
        {
            const auto& next_N_deme = current_deme_parameters.next_deme_sizes.get();
            for (std::size_t i = 0; i < next_N_deme.size(); ++i)
                {
                    if (next_N_deme[i] > 0 && fitnesses.lookups[i] == nullptr)
                        {
                            if (M.empty())
                                {
                                    detail::no_valid_parents(i, generation,
                                                             next_N_deme[i]);
                                }
                            else
                                {
                                    detail::check_migration_in(i, generation, M);
                                }
                        }
                }
        }

    } // namespace discrete_demography
} // namespace fwdpy11

#endif
