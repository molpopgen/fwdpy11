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

#ifndef FWDPY11_DEME_PROPERTY_TYPES_HPP
#define FWDPY11_DEME_PROPERTY_TYPES_HPP

#include <cstdint>
#include <vector>
#include "detail.hpp"
#include <fwdpp/named_type.hpp>

namespace fwdpy11
{
    namespace discrete_demography
    {
        using current_deme_sizes_vector = fwdpp::strong_types::named_type<
            std::vector<std::uint32_t>, detail::current_deme_sizes_vector_t>;
        using next_deme_sizes_vector = fwdpp::strong_types::named_type<
            std::vector<std::uint32_t>, detail::next_deme_sizes_vector_t>;
        using growth_rates_onset_times_vector
            = fwdpp::strong_types::named_type<
                std::vector<std::uint32_t>,
                detail::growth_onset_times_vector_t>;
        using growth_initial_size_vector = fwdpp::strong_types::named_type<
            std::vector<std::uint32_t>, detail::growth_initial_size_vector_t>;
        using growth_rates_vector
            = fwdpp::strong_types::named_type<std::vector<double>,
                                              detail::growth_rates_vector_t>;
        using selfing_rates_vector
            = fwdpp::strong_types::named_type<std::vector<double>,
                                              detail::selfing_rates_vector_t>;
    } // namespace discrete_demography
} // namespace fwdpy11

#endif

