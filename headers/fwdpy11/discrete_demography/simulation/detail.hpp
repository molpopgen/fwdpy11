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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_DETAIL_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_DETAIL_HPP

#include <utility>

namespace fwdpy11
{
    namespace discrete_demography
    {
        namespace detail
        {
            using event_range = std::pair<std::size_t, const std::size_t>;

            struct mass_migration_range_t
            {
            };
            struct deme_size_change_range_t
            {
            };
            struct growth_rate_change_range_t
            {
            };
            struct selfing_rate_change_range_t
            {
            };
            struct migration_rate_change_range_t
            {
            };

            struct current_deme_sizes_vector_t
            {
            };
            struct next_deme_sizes_vector_t
            {
            };
            struct growth_onset_times_vector_t
            {
            };
            struct growth_initial_size_vector_t
            {
            };
            struct growth_rates_vector_t
            {
            };
            struct selfing_rates_vector_t
            {
            };
        } // namespace detail
    }     // namespace discrete_demography
} // namespace fwdpy11

#endif
