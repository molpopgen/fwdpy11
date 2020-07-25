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
#include <fwdpy11/genetic_value_data/genetic_value_data.hpp>
#include <pybind11/pybind11.h>

using genetic_values_buffer_proxy = fwdpy11::double_array_proxy;

struct GeneticValueIsTraitData
{
    pybind11::list
    fill_list(pybind11::object o1, pybind11::object o2)
    {
        pybind11::list rv;
        rv.append(o1);
        rv.append(o2);
        return rv;
    }
    
    fwdpy11::DiploidMetadata offspring_metadata_copy, parent1_copy, parent2_copy;
    pybind11::object offspring_metadata, parent1, parent2;
    genetic_values_buffer_proxy buffer;
    pybind11::object genetic_values;
    pybind11::list parental_metadata;
    std::size_t offspring_metadata_index;

    GeneticValueIsTraitData()
        : offspring_metadata_copy{}, parent1_copy{}, parent2_copy{},
          offspring_metadata{
              pybind11::cast<fwdpy11::DiploidMetadata*>(&offspring_metadata_copy)},
          parent1{pybind11::cast<fwdpy11::DiploidMetadata*>(&parent1_copy)},
          parent2{pybind11::cast<fwdpy11::DiploidMetadata*>(&parent2_copy)}, buffer{},
          genetic_values{pybind11::cast<genetic_values_buffer_proxy*>(&buffer)},
          parental_metadata{fill_list(parent1, parent2)},
          offspring_metadata_index{std::numeric_limits<std::size_t>::max()}
    {
    }
};

inline void
set_data(const fwdpy11::DiploidGeneticValueToFitnessData& input_data,
         GeneticValueIsTraitData& data)
{
    data.offspring_metadata_copy = input_data.offspring_metadata.get();
    data.buffer.data = const_cast<double*>(input_data.gvalues.get().data());
    data.buffer.size = input_data.gvalues.get().size();
    data.offspring_metadata_index = input_data.metadata_index;
    data.parent1_copy = input_data.parent1_metadata.get();
    data.parent2_copy = input_data.parent2_metadata.get();
}

#endif
