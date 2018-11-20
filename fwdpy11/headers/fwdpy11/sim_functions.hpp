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
#ifndef FWDPY11_SIM_FUNCTIONS_HPP_
#define FWDPY11_SIM_FUNCTIONS_HPP_

/*!
  \file sim_features.hpp

  fwdpy11 serves as a test-bed for fwdpp.  Sometimes, I realize that fwdpp
  is missing a feature.  Those features may first appear here before getting
  moved over to fwdpp.
*/

#include <tuple>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <fwdpp/util.hpp>

namespace fwdpy11
{
    /// This function is similar in name and interface to the current fwdpp
    /// function in fwdpp/util.hpp. However, it accepts an additional bool
    /// to determine whether or not fixed, non-neutral mutations are flagged
    /// for recycling or not.
    ///
    /// It differs from the current fwdpp version in that:
    /// 1. It uses std::lower_bound to make sure that fixations/fixation times
    /// are sorted by position
    /// 2. It uses binary searches (again, lower_bound) to guard against
    /// re-inserting the same non-neutral fixation over and over.
    ///
    /// The reason for these changes is that the use case is sims of
    /// phenotypes.
    /// We keep fixations in the pop so that they contribute to trait values.
    /// Thus, w/o the
    /// searches, we'd keep re-inserting a fixation each generation.
    ///
    /// \note: lookup must be compatible with
    /// lookup->erase(lookup->find(double))
    template <typename mcont_t, typename fixation_container_t,
              typename fixation_time_container_t,
              typename mutation_lookup_table>
    void
    update_mutations(mcont_t &mutations, fixation_container_t &fixations,
                     fixation_time_container_t &fixation_times,
                     mutation_lookup_table &lookup,
                     std::vector<fwdpp::uint_t> &mcounts,
                     const unsigned &generation, const unsigned &twoN,
                     const bool remove_selected_fixations)
    {
        using namespace fwdpp;
        static_assert(
            typename traits::is_mutation_t<
                typename mcont_t::value_type>::type(),
            "mutation_type must be derived from fwdpp::mutation_base");
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error("DEBUG: container size error");
            }
#endif
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
#ifndef NDEBUG
                if (mcounts[i] > twoN)
                    {
                        throw std::runtime_error(
                            "DEBUG: mutation count too large");
                    }
#endif
                if (mcounts[i] == twoN)
                    {
                        auto loc = std::lower_bound(
                            fixations.begin(), fixations.end(),
                            std::make_tuple(mutations[i].pos, mutations[i].g),
                            [](const typename fixation_container_t::value_type
                                   &mut,
                               const std::tuple<double, std::uint32_t>
                                   &value) noexcept {
                                return std::tie(mut.pos, mut.g) < value;
                            });
                        auto d = std::distance(fixations.begin(), loc);
                        if (mutations[i].neutral
                            || remove_selected_fixations == true)
                            {
                                fixations.insert(loc, mutations[i]);
                                fixation_times.insert(
                                    fixation_times.begin() + d, generation);
                                mcounts[i] = 0; // set count to zero to mark
                                                // mutation as "recyclable"
                                auto itr
                                    = lookup.equal_range(mutations[i].pos);
                                // Make position max double so that a user
                                // cannot accidentally track this as a zero-frequency
                                // variant
                                mutations[i].pos
                                    = std::numeric_limits<double>::max();
                                while (itr.first != itr.second)
                                    {
                                        if (itr.first->second == i)
                                            {
                                                lookup.erase(itr.first);
                                                break;
                                            }
                                        ++itr.first;
                                    }
                            }
                        else
                            {
                                if (loc == fixations.end()
                                    || (loc->pos != mutations[i].pos
                                        && loc->g != mutations[i].g))
                                    {
                                        fixations.insert(loc, mutations[i]);
                                        fixation_times.insert(
                                            fixation_times.begin() + d,
                                            generation);
                                    }
                            }
                    }
                else if (!mcounts[i])
                    {
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        // Make position max double so that a user
                                        // cannot accidentally track this as a zero-frequency
                                        // variant
                                        mutations[itr.first->second].pos
                                            = std::numeric_limits<
                                                double>::max();
                                        lookup.erase(itr.first);
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
            }

        //    struct update_mutations_wrapper
        //    {
        //        template <typename... args>
        //        inline void
        //        operator()(args &&... Args) const
        //        {
        //            fwdpp::update_mutations(std::forward<args>(Args)...);
        //        }
        //    };
        //
        //    struct update_mutations_n_wrapper
        //    {
        //        template <typename... args>
        //        inline void
        //        operator()(args &&... Args) const
        //        {
        //            fwdpy11::update_mutations_n(std::forward<args>(Args)...);
        //        }
        //    };
    }
} // namespace fwdpy11

#endif
