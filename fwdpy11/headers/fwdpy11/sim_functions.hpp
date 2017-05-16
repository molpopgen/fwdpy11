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
#include <algorithm>
#include <fwdpp/util.hpp>

namespace fwdpy11
{
    /*!
          Label all fixed neutral variant and all extinct variants for
       recycling.
          Copy fixations and fixation times
          for neutral mutations into containers.

          Fixed, non-neutral variants get copied into fixations and fixation
       times,
          so that fixation times
          can get records.

          This function is identical in name and interface to the current fwdpp
          function in fwdpp/util.hpp.

          It differs from the current fwdpp version in that:
          1. It uses std::lower_bound to make sure that fixations/fixation
       times
          are sorted by position
          2. It uses binary searches (again, lower_bound) to guard against
          re-inserting the same
          non-neutral fixation over and over.

          The reason for these changes is that the use case is sims of
       phenotypes.
          We keep fixations
          in the pop so that they contribute to trait values.  Thus, w/o the
          searches, we'd keep re-inserting
          a fixation each generation.

          \note: lookup must be compatible with
       lookup->erase(lookup->find(double))
        */
    template <typename mcont_t, typename fixation_container_t,
              typename fixation_time_container_t,
              typename mutation_lookup_table>
    void
    update_mutations_n(mcont_t &mutations, fixation_container_t &fixations,
                       fixation_time_container_t &fixation_times,
                       mutation_lookup_table &lookup,
                       std::vector<KTfwd::uint_t> &mcounts,
                       const unsigned &generation, const unsigned &twoN)
    {
        using namespace KTfwd;
        static_assert(
            typename traits::is_mutation_t<
                typename mcont_t::value_type>::type(),
            "mutation_type must be derived from KTfwd::mutation_base");
        assert(mcounts.size() == mutations.size());
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
                assert(mcounts[i] <= twoN);
                if (mcounts[i] == twoN)
                    {
                        auto loc = std::lower_bound(
                            fixations.begin(), fixations.end(),
                            mutations[i].pos,
                            [](const typename fixation_container_t::value_type
                                   &__mut,
                               const double &__value) noexcept {
                                return __mut.pos < __value;
                            });
                        auto d = std::distance(fixations.begin(), loc);
                        if (mutations[i].neutral)
                            {
                                fixations.insert(loc, mutations[i]);
                                fixation_times.insert(
                                    fixation_times.begin() + d, generation);
                                mcounts[i] = 0; // set count to zero to mark
                                                // mutation as "recyclable"
                                lookup.erase(mutations[i].pos); // remove
                                                                // mutation
                                                                // position
                                                                // from lookup
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
                if (!mcounts[i])
                    lookup.erase(mutations[i].pos);
            }
    }

    struct update_mutations_wrapper
    {
        template <typename... args>
        inline void
        operator()(args &&... Args) const
        {
            KTfwd::update_mutations(std::forward<args>(Args)...);
        }
    };

    struct update_mutations_n_wrapper
    {
        template <typename... args>
        inline void
        operator()(args &&... Args) const
        {
            fwdpy11::update_mutations_n(std::forward<args>(Args)...);
        }
    };
}

#endif
