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
#include "../MassMigration.hpp"
#include "../DiscreteDemography.hpp"
#include "deme_properties.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename METADATATYPE>
        inline void
        mass_migration(
            const GSLrng_t& rng, std::uint32_t t,
            mass_migration_range& mass_migration_ranges,
            growth_rates_vector& growth_rates,
            growth_rates_onset_times_vector& growth_rate_onset_times,
            growth_initial_size_vector& growth_initial_sizes,
            std::vector<METADATATYPE>& metadata)
        {
            auto& ref = mass_migration_ranges.get();
            if (ref.first < ref.second && t < ref.first->when)
                {
                    return;
                }
            while (ref.first < ref.second && ref.first->when == t)
                {
                    ref.first = apply_mass_migrations(
                        rng, t, ref.first, ref.second, growth_rates,
                        growth_rate_onset_times, growth_initial_sizes,
                        metadata);
                }
        }

        namespace detail
        {
            inline void
            update_current_deme_sizes(const std::uint32_t t,
                                      deme_size_change_range& event_ranges,
                                      deme_properties& sizes_rates)
            // NOTE: this function may resent growth rates to zero.
            // see SetDemeSize for details.
            {
                deme_size_change_range::value_type& range = event_ranges.get();
                if (range.first < range.second && t < range.first->when)
                    {
                        return;
                    }
                current_deme_sizes_vector::value_type& current_deme_sizes
                    = sizes_rates.current_deme_sizes.get();
                growth_rates_vector::value_type& growth_rates
                    = sizes_rates.growth_rates.get();
                growth_rates_onset_times_vector::value_type&
                    growth_rate_onset_times
                    = sizes_rates.growth_rate_onset_times.get();
                growth_initial_size_vector::value_type& growth_initial_sizes
                    = sizes_rates.growth_initial_sizes.get();
                for (; range.first < range.second && range.first->when == t;
                     ++range.first)
                    {
                        auto deme = range.first->deme;
                        current_deme_sizes[deme] = range.first->new_size;
                        if (range.first->resets_growth_rate == true)
                            {
                                growth_rates[deme] = NOGROWTH;
                            }
                        // Deme size has changed, so we reset
                        // the onset of growth and the initial N for this deme
                        // NOTE: this needs careful documentation!
                        growth_rate_onset_times[deme] = t;
                        growth_initial_sizes[deme] = current_deme_sizes[deme];
                    }
            }

            inline void
            update_growth_rates(
                const std::uint32_t t,
                growth_rate_change_range& growth_rate_change_tracker,
                deme_properties& sizes_rates)
            {
                auto& range = growth_rate_change_tracker.get();
                if (range.first < range.second && t < range.first->when)
                    {
                        return;
                    }
                auto& rates = sizes_rates.growth_rates.get();
                auto& onsets = sizes_rates.growth_rate_onset_times.get();
                auto& current_deme_sizes
                    = sizes_rates.current_deme_sizes.get();
                auto& N0 = sizes_rates.growth_initial_sizes.get();
                auto& Ncurr = sizes_rates.current_deme_sizes.get();
                for (; range.first < range.second && range.first->when == t;
                     ++range.first)
                    {
                        auto deme = range.first->deme;
                        if (range.first->G != NOGROWTH && Ncurr[deme] == 0)
                            {
                                throw DemographyError(
                                    "attempt to change growth rate in extinct "
                                    "deme");
                            }
                        rates[deme] = range.first->G;
                        onsets[deme] = t;
                        N0[deme] = current_deme_sizes[deme];
                    }
            }

            inline void
            update_selfing_rates(
                const std::uint32_t t,
                selfing_rate_change_range& selfing_rate_change_tracker,
                deme_properties& sizes_rates)
            {
                auto& range = selfing_rate_change_tracker.get();
                if (range.first < range.second && t < range.first->when)
                    {
                        return;
                    }
                auto& rates = sizes_rates.selfing_rates.get();
                auto& Ncurr = sizes_rates.current_deme_sizes.get();
                for (; range.first < range.second && range.first->when == t;
                     ++range.first)
                    {
                        if (Ncurr[range.first->deme] == 0)
                            {
                                throw DemographyError("attempt to set selfing "
                                                      "rate in extinct deme");
                            }
                        rates[range.first->deme] = range.first->S;
                    }
            }

            inline void
            update_migration_matrix(
                const std::uint32_t t,
                migration_rate_change_range& migration_rate_change_tracker,
                std::unique_ptr<MigrationMatrix>& M)

            {
                auto& range = migration_rate_change_tracker.get();
                if (range.first < range.second && t < range.first->when)
                    {
                        return;
                    }
                for (; range.first < range.second && range.first->when == t;
                     ++range.first)
                    {
                        M->set_migration_rates(range.first->deme,
                                               range.first->migrates);
                    }
            }

            inline std::uint32_t
            apply_growth_rates_get_next_global_N(const std::uint32_t t,
                                                 deme_properties& sizes_rates)
            {
                auto& Ncurr = sizes_rates.current_deme_sizes.get();
                auto& Nnext = sizes_rates.next_deme_sizes.get();
                auto& G = sizes_rates.growth_rates.get();
                auto& N0 = sizes_rates.growth_initial_sizes.get();
                auto& onset = sizes_rates.growth_rate_onset_times.get();
                std::uint32_t next_global_N = 0;
                for (std::size_t deme = 0; deme < Nnext.size(); ++deme)
                    {
                        if (G[deme] != NOGROWTH)
                            {
                                if (Ncurr[deme] == 0)
                                    {
                                        throw DemographyError(
                                            "growth is happening in an "
                                            "extinct deme");
                                    }
                                double next_size = std::round(
                                    static_cast<double>(N0[deme])
                                    * std::pow(G[deme],
                                               static_cast<double>(
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
        } // namespace detail

        inline std::uint32_t
        apply_demographic_events(std::uint32_t t,
                                 DiscreteDemography& demography,
                                 std::unique_ptr<MigrationMatrix>& M,
                                 deme_properties& sizes_rates)
        {
            std::copy(begin(sizes_rates.current_deme_sizes.get()),
                      end(sizes_rates.current_deme_sizes.get()),
                      begin(sizes_rates.next_deme_sizes.get()));
            // Step 1, do the discrete changes of deme sizes
            // NOTE: this may reset growth rates to zero
            detail::update_current_deme_sizes(
                t, demography.deme_size_change_tracker, sizes_rates);

            // Step 2: set new growth rates
            detail::update_growth_rates(
                t, demography.growth_rate_change_tracker, sizes_rates);
            // Step 3: update selfing rates
            detail::update_selfing_rates(
                t, demography.selfing_rate_change_tracker, sizes_rates);
            detail::update_migration_matrix(
                t, demography.migration_rate_change_tracker, M);
            // Step 4: set next deme sizes and apply growth rates
            std::copy(begin(sizes_rates.current_deme_sizes.get()),
                      end(sizes_rates.current_deme_sizes.get()),
                      begin(sizes_rates.next_deme_sizes.get()));
            auto next_global_N
                = detail::apply_growth_rates_get_next_global_N(t, sizes_rates);
            return next_global_N;
        }

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
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
