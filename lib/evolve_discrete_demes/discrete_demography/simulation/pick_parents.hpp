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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_PICK_PARENTS_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_PICK_PARENTS_HPP

#include <cstdint>
#include <sstream>
#include <fwdpy11/rng.hpp>
#include "core/demes/forward_graph.hpp"
#include "fwdpp/gsl_discrete.hpp"
#include "multideme_fitness_lookups.hpp"
#include "mating_event_type.hpp"

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        struct parent_data
        {
            std::size_t parent1, parent2;
            std::int32_t deme1, deme2;
            mating_event_type mating;

            parent_data(std::size_t p1, std::size_t p2, std::int32_t d1, std::int32_t d2,
                        mating_event_type m)
                : parent1(p1), parent2(p2), deme1(d1), deme2(d2), mating(m)
            {
            }
        };

        parent_data
        pick_parents(const fwdpy11::GSLrng_t& rng, const std::int32_t offspring_deme,
                     const fwdpy11_core::ForwardDemesGraph& demography,
                     const fwdpp::gsl_ran_discrete_t_ptr& ancestor_deme_lookup,
                     const multideme_fitness_bookmark& fitness_bookmark,
                     const multideme_fitness_lookups<std::uint32_t>& wlookups,
                     const bool allow_residual_selfing);

    } // namespace discrete_demography
} // namespace fwdpy11_core

#endif
