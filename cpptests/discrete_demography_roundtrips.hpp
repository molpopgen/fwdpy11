#pragma once
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
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/discrete_demography/simulation/mating_event_type.hpp>

struct MatingEventRecord
{
    std::uint32_t generation;
    std::int32_t parental_deme, offspring_deme;
    fwdpy11::discrete_demography::mating_event_type mating_event;

    MatingEventRecord(std::uint32_t g, std::int32_t pd, std::int32_t od,
                      fwdpy11::discrete_demography::mating_event_type mt);
};

std::unordered_map<int, int>
get_deme_sizes(const std::vector<fwdpy11::DiploidMetadata>&);

std::vector<MatingEventRecord> DiscreteDemography_roundtrip(
    const fwdpy11::GSLrng_t& rng, fwdpy11::DiploidPopulation& pop,
    fwdpy11::discrete_demography::DiscreteDemography& demography, std::uint32_t ngens);
