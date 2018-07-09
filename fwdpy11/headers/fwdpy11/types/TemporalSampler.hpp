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
#ifndef FWDPY11_SAMPLERS_HPP__
#define FWDPY11_SAMPLERS_HPP__

#include <functional>
#include <string>
#include <memory>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

namespace fwdpy11
{
    // The following two typedefs are used for the case
    // of bare C++ samplers, or pure Python functions or 
    // callable classes
    using SlocusPop_temporal_sampler
        = std::function<void(const fwdpy11::SlocusPop&)>;

    using MlocusPop_temporal_sampler
        = std::function<void(const fwdpy11::MlocusPop&)>;

    struct TemporalSampler
    {
        virtual std::string serialize() const = 0;
        virtual void deserialize(std::string) = 0;
        virtual std::unique_ptr<TemporalSampler> clone() const = 0;
    };
}

#endif
