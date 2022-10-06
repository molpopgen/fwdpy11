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
#pragma once

#include <vector>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/rng.hpp>

struct new_mutation_data
{
    double effect_size, dominance;
    std::vector<double> esizes, heffects;
    decltype(fwdpy11::Mutation::xtra) label;

    new_mutation_data(double e, double h, std::vector<double> esizes,
                      std::vector<double> heffects, decltype(fwdpy11::Mutation::xtra) l);
};

std::size_t add_mutation(const fwdpy11::GSLrng_t& rng, const double left,
                         const double right, const fwdpp::ts::table_index_t ndescendants,
                         const fwdpp::ts::table_index_t deme,
                         const new_mutation_data& data, fwdpy11::DiploidPopulation& pop);
