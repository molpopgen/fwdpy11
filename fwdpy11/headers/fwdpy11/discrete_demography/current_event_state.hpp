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

#ifndef FWDPY11_DISCRETE_CURRENT_EVENT_STATE_HPP
#define FWDPY11_DISCRETE_CURRENT_EVENT_STATE_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <utility>

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename T> struct current_event_state
        {
            std::vector<T> events;
            std::pair<std::size_t, std::size_t> event_range;
            template <typename Input>
            current_event_state(Input&& input)
                : events(std::forward<Input>(input)), event_range{0, events.size()}
            {
                if (!std::is_sorted(begin(events), end(events)))
                    {
                        throw std::invalid_argument("events list is not sorted");
                    }
            }

            const T&
            event(std::size_t index) const
            {
                return events[index];
            }

            const T&
            event() const
            {
                return event(current());
            }

            auto
            when() const
            {
                return event().when;
            }

            std::size_t
            current() const
            {
                return event_range.first;
            }

            std::size_t&
            current()
            {
                return event_range.first;
            }

            std::size_t
            last() const
            {
                return event_range.second;
            }
        };
    }
}

#endif
