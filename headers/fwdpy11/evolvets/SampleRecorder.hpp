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
#ifndef FWDPY11_EVOLVETS_SAMPLE_RECORDER_HPP
#define FWDPY11_EVOLVETS_SAMPLE_RECORDER_HPP

#include <stdexcept>
#include <vector>
#include <fwdpp/forward_types.hpp>

namespace fwdpy11
{
    struct SampleRecorder
    /*! This type will be passed to a user-supplied callable
     * during simulations with tree sequences.
     *
     * The callable is expected to populate the object with 
     * the indexes of individuals (NOT NODES) to be preserved
     * as ancient samples.
     *
     * It is then the responsibility of the simulation routine
     * to map individual indexes to node indexes and preserve
     * them as ancient samples, along with relevant metadata, etc.
     */
    {
        std::vector<fwdpp::uint_t> samples;

        SampleRecorder() : samples{}
        {
        }

        void
        add_sample(const fwdpp::uint_t i)
        {
            samples.push_back(i);
        }
    };
} // namespace fwdpy11
#endif
