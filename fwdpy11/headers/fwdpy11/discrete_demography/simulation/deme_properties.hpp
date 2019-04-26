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

#ifndef FWDPY11_DEME_PROPERTIES_HPP
#define FWDPY11_DEME_PROPERTIES_HPP

#include "../constants.hpp"
#include "deme_property_types.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct deme_properties
        {
            template <typename METADATATYPE>
            current_deme_sizes_vector
            init_deme_sizes_from_metadata(
                std::uint32_t max_number_demes,
                const std::vector<METADATATYPE>& metadata)
            {
                current_deme_sizes_vector v(
                    current_deme_sizes_vector::value_type(max_number_demes,
                                                          0));
                auto& ref = v.get();
                std::fill(begin(ref), end(ref), 0);
                for (auto&& i : metadata)
                    {
                        ref[i.deme]++;
                    }
                return v;
            }

            current_deme_sizes_vector current_deme_sizes;
            next_deme_sizes_vector next_deme_sizes;
            growth_rates_onset_times_vector growth_rate_onset_times;
            growth_initial_size_vector growth_initial_sizes;
            growth_rates_vector growth_rates;
            selfing_rates_vector selfing_rates;

            template <typename METADATATYPE>
            deme_properties(std::uint32_t max_number_demes,
                            const std::vector<METADATATYPE>& metadata)
                : current_deme_sizes(
                    init_deme_sizes_from_metadata(max_number_demes, metadata)),
                  next_deme_sizes(
                      next_deme_sizes_vector::value_type(max_number_demes, 0)),
                  growth_rate_onset_times(
                      growth_rates_onset_times_vector::value_type(
                          max_number_demes, 0)),
                  growth_initial_sizes(current_deme_sizes.get()),
                  growth_rates(growth_rates_vector::value_type(
                      max_number_demes, NOGROWTH)),
                  selfing_rates(
                      selfing_rates_vector::value_type(max_number_demes, 0.))
            {
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
