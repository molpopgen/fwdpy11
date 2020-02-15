//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
#ifndef FWDPY11_SET_MODEL_STATE_HPP
#define FWDPY11_SET_MODEL_STATE_HPP

#include <vector>
#include <memory>
#include "../DiscreteDemography.hpp"
#include "demographic_model_state.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {

        template <typename METADATATYPE>
        inline demographic_model_state_pointer
        initialize_model_state(std::uint32_t generation,
                               const std::vector<METADATATYPE>& metadata,
                               DiscreteDemography& demography)
        {
            // "steal" pointer from input
            auto rv = demography.get_model_state();

            if (rv == nullptr || generation == 0)
            // If there is no state, then we need to make
            // one.  If there is a state, but the generation
            // is zero, then we assume that the demography
            // has been used for a different simulatin replicate
            // and thus reset it.
                {
                    demography.update_event_times(generation);
                    rv.reset(
                        new demographic_model_state(metadata, demography));
                }
            return rv;
        }

        inline void
        save_model_state(demographic_model_state_pointer& state,
                         DiscreteDemography& demography)
        {
            demography.set_model_state(state);
        }
    } // namespace discrete_demography
} // namespace fwdpy11
#endif
