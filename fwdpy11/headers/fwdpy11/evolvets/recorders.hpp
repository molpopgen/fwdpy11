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
#ifndef FWDPY11_EVOLVETS_RECORDERS_HPP
#define FWDPY11_EVOLVETS_RECORDERS_HPP

#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include "samplerecorder.hpp"

namespace fwdpy11
{
    struct no_ancient_samples
    /*! When no ancient samples are tracked, 
     * this will provide the most efficient 
     * way to "do nothing" b/c the 
     * repeated C++/Py round trips are avoided.
     */
    {
        inline void
        operator()(const fwdpy11::SlocusPop&,
                   const fwdpy11::samplerecorder&) const
        {
        }
        inline void
        operator()(const fwdpy11::MlocusPop&,
                   const fwdpy11::samplerecorder&) const
        {
        }
    };
} // namespace fwdpy11
#endif
