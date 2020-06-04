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

#ifndef FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_MANAGER_HPP
#define FWDPY11_DISCRETE_DEMOGRAPHY_DISCRETE_DEMOGRAPHY_MANAGER_HPP

#include <memory>
#include <cstdint>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include "get_max_number_of_demes.hpp"
#include "multideme_fitness_lookups.hpp"
#include "migration_lookup.hpp"
#include "deme_properties.hpp"
#include "../DiscreteDemography.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class demographic_model_state
        /// Added in 0.6.0 to hold and manage
        /// the relevant data structures.
        {
          private:
            std::unique_ptr<MigrationMatrix>
            init_migmatrix(const std::unique_ptr<const MigrationMatrix>& Minput)
            {
                if (Minput == nullptr)
                    {
                        return nullptr;
                    }
                return std::unique_ptr<MigrationMatrix>(new MigrationMatrix(*Minput));
            }

            std::uint32_t next_global_N;

          public:
            const std::int32_t maxdemes;
            multideme_fitness_lookups<std::uint32_t> fitnesses;
            deme_properties sizes_rates;
            std::unique_ptr<MigrationMatrix> M;
            migration_lookup miglookup;

            // NOTE: demography.update_event_times() needs to have been
            // called first!
            template <typename METADATATYPE>
            demographic_model_state(const std::vector<METADATATYPE>& metadata,
                                    DiscreteDemography& demography)
                : next_global_N(0),
                  maxdemes(get_max_number_of_demes()(metadata, demography)),
                  fitnesses(maxdemes), sizes_rates(maxdemes, metadata),
                  M(init_migmatrix(demography.migmatrix)),
                  miglookup(maxdemes, M == nullptr)
            {
            }

            // This constructor is only used when resetting
            // the state from an event like pickling a DiscreteDemography
            // instance.
            demographic_model_state(std::int32_t maxdemes_, deme_properties sizes_rates_,
                                    std::unique_ptr<MigrationMatrix> M_)
                : next_global_N(0), maxdemes(maxdemes_), fitnesses(maxdemes),
                  sizes_rates(std::move(sizes_rates_)), M(std::move(M_)),
                  miglookup(maxdemes, M == nullptr)
            {
                next_global_N
                    = std::accumulate(begin(sizes_rates.next_deme_sizes.get()),
                                      end(sizes_rates.next_deme_sizes.get()), 0u);
            }

            void
            set_next_global_N(std::uint32_t N)
            {
                next_global_N = N;
            }

            std::uint32_t
            ttlN_next() const
            {
                return next_global_N;
            }

            bool
            will_go_globally_extinct() const
            {
                return ttlN_next() == 0;
            }
        };

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
                    rv.reset(new demographic_model_state(metadata, demography));
                }
            return rv;
        }

        inline void
        save_model_state(demographic_model_state_pointer state,
                         DiscreteDemography& demography)
        {
            demography.set_model_state(std::move(state));
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
