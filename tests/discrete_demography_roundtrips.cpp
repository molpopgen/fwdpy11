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
#include <utility>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <tuple>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/discrete_demography/simulation.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::DiploidMetadata>);

namespace ddemog = fwdpy11::discrete_demography;

// The return value records migration events: offspring generation,
// parent deme, offspring deme.
using roundtrip_rv_type
    = std::vector<std::tuple<std::uint32_t, std::int32_t, std::int32_t,
                             ddemog::mating_event_type>>;

roundtrip_rv_type
DiscreteDemography_roundtrip(const fwdpy11::GSLrng_t& rng,
                             fwdpy11::DiploidPopulation& pop,
                             ddemog::DiscreteDemography& demography,
                             std::uint32_t ngens)
// Assumptions:
// 1. Initial deme labels are set by the user. NOTE: validated by manager object
{
    demography.update_event_times(pop.generation);
    ddemog::discrete_demography_manager ddemog_manager(pop.diploid_metadata,
                                                       demography);
    decltype(pop.diploid_metadata) offspring_metadata;
    offspring_metadata.reserve(pop.N);
    roundtrip_rv_type rv;
    for (std::uint32_t gen = 0; gen < ngens; ++gen, ++pop.generation)
        {
            ddemog::update_demography_manager(rng, pop.generation,
                                              pop.diploid_metadata, demography,
                                              ddemog_manager);

            // Generate the offspring
            for (std::int32_t deme = 0; deme < ddemog_manager.maxdemes; ++deme)
                {
                    for (unsigned ind = 0;
                         ind < ddemog_manager.sizes_rates.next_deme_sizes
                                   .get()[deme];
                         ++ind)
                        {
                            auto pdata = ddemog::pick_parents(
                                rng, deme, ddemog_manager.miglookup,
                                ddemog_manager.sizes_rates.current_deme_sizes,
                                ddemog_manager.sizes_rates.selfing_rates,
                                ddemog_manager.fitnesses);
                            if (pdata.deme1 != deme)
                                {
                                    rv.emplace_back(pop.generation + 1,
                                                    pdata.deme1, deme,
                                                    pdata.mating);
                                }
                            if (pdata.deme2 != deme)
                                {
                                    rv.emplace_back(pop.generation + 1,
                                                    pdata.deme2, deme,
                                                    pdata.mating);
                                }
                            offspring_metadata.emplace_back(
                                fwdpy11::DiploidMetadata{
                                    0.0,
                                    0.,
                                    1.,
                                    { 0, 0, 0 },
                                    offspring_metadata.size(),
                                    { pdata.parent1, pdata.parent2 },
                                    deme,
                                    -1,
                                    { -1, -1 } });
                        }
                }
            pop.diploid_metadata.swap(offspring_metadata);
            offspring_metadata.clear();

            pop.N = static_cast<std::uint32_t>(pop.diploid_metadata.size());
        }
    return rv;
}

PYBIND11_MODULE(discrete_demography_roundtrips, m)
{
    m.def("DiscreteDemography_roundtrip", &DiscreteDemography_roundtrip);

    pybind11::enum_<ddemog::mating_event_type>(m, "MatingEventType",
                                               pybind11::arithmetic())
        .value("outcrossing", ddemog::mating_event_type::outcrossing)
        .value("selfing", ddemog::mating_event_type::selfing)
        .value("cloning", ddemog::mating_event_type::cloning);
}
