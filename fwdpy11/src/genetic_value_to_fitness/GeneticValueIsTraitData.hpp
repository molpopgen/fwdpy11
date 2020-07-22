//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWPY11_GENETIC_VALUE_IS_TRAIT_DATA_HPP
#define FWPY11_GENETIC_VALUE_IS_TRAIT_DATA_HPP

#include <limits>
#include <fwdpy11/util/array_proxy.hpp>
#include <fwdpy11/types/Diploid.hpp>
#include <pybind11/pybind11.h>

using genetic_values_buffer_proxy = fwdpy11::double_array_proxy;

struct GeneticValueIsTraitData
{
    fwdpy11::DiploidMetadata offspring_metadata_copy;
    pybind11::object offspring_metadata;
    genetic_values_buffer_proxy buffer;
    pybind11::object genetic_values;
    std::size_t offspring_metadata_index;

    GeneticValueIsTraitData()
        : offspring_metadata_copy{},
          offspring_metadata{
              pybind11::cast<fwdpy11::DiploidMetadata*>(&offspring_metadata_copy)},
          buffer{}, genetic_values{pybind11::cast<genetic_values_buffer_proxy*>(
                        &buffer)},
          offspring_metadata_index{std::numeric_limits<std::size_t>::max()}
    {
    }
};

#endif
