//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
#ifndef FWDPY11_UPDATE_DEMOGRAPHY_MANAGER_HPP
#define FWDPY11_UPDATE_DEMOGRAPHY_MANAGER_HPP

#include <memory>
#include <fwdpy11/rng.hpp>
#include "discrete_demography_manager.hpp"
#include "functions.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename METADATATYPE>
        inline void
        update_demography_manager(
            const GSLrng_t &rng, const std::uint32_t generation,
            std::vector<METADATATYPE> &metadata,
            DiscreteDemography &demography,
            std::unique_ptr<discrete_demography_manager> &ddemog_manager)
        {
            mass_migration(rng, generation, demography.mass_migration_tracker,
                           ddemog_manager->sizes_rates.growth_rates,
                           ddemog_manager->sizes_rates.growth_rate_onset_times,
                           ddemog_manager->sizes_rates.growth_initial_sizes,
                           metadata);
            get_current_deme_sizes(
                metadata, ddemog_manager->sizes_rates.current_deme_sizes);
            ddemog_manager->fitnesses.update(
                ddemog_manager->sizes_rates.current_deme_sizes, metadata);
            auto next_global_N = apply_demographic_events(
                generation, demography, ddemog_manager->M,
                ddemog_manager->sizes_rates);
            ddemog_manager->set_next_global_N(next_global_N);
            build_migration_lookup(
                ddemog_manager->M,
                ddemog_manager->sizes_rates.current_deme_sizes,
                ddemog_manager->sizes_rates.selfing_rates,
                ddemog_manager->miglookup);
        }
    } // namespace discrete_demography
} // namespace fwdpy11
#endif
