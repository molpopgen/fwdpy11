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

#ifndef FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_FUNCTIONS_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_SIMULATION_FUNCTIONS_HPP

#include <vector>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <core/demes/forward_graph.hpp>
#include "multideme_fitness_lookups.hpp"

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        void validate_parental_state(
            std::uint32_t generation,
            const multideme_fitness_lookups<std::uint32_t>& fitnesses,
            const fwdpy11_core::ForwardDemesGraph& demography);

    } // namespace discrete_demography
} // namespace fwdpy11_core

#endif
