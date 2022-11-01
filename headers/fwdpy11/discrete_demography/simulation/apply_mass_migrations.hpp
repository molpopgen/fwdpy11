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
#ifndef FWDPY11_APPLY_MOVE_OR_COPY_HPP
#define FWDPY11_APPLY_MOVE_OR_COPY_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <stack>
#include <unordered_map>
#include <sstream>
#include <gsl/gsl_randist.h>
#include "../../rng.hpp"
#include "../MassMigration.hpp"
#include "../constants.hpp"
#include "../exceptions.hpp"
#include "../current_event_state.hpp"
#include "deme_property_types.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        namespace detail
        {
            using deme_map_t
                = std::unordered_map<std::int32_t, std::vector<std::size_t>>;
            template <typename METADATATYPE>
            inline deme_map_t
            build_deme_map(const std::vector<METADATATYPE>& metadata)
            {
                deme_map_t rv;
                for (std::size_t i = 0; i < metadata.size(); ++i)
                    {
                        if (i != metadata[i].label)
                            {
                                throw std::runtime_error(
                                    "metadata label does not equal index");
                            }
                        rv[metadata[i].deme].push_back(i);
                    }
                return rv;
            }

            //using move_stack = std::stack<std::size_t, std::vector<std::size_t>>;
            using vector_backed_move_stack
                = std::stack<std::size_t, std::vector<std::size_t>>;

            struct move_stack
            {
                std::size_t initial_N;
                vector_backed_move_stack ms;
                template <typename T>
                move_stack(std::size_t N, T&& t) : initial_N(N), ms(std::forward<T>(t))
                {
                }
            };

            using move_map_t = std::unordered_map<std::int32_t, move_stack>;

            inline move_map_t
            build_move_sources(const GSLrng_t& rng, const deme_map_t& deme_map)
            {
                move_map_t rv;

                for (auto&& m : deme_map)
                    {
                        auto indexes(m.second);
                        if (indexes.empty())
                            {
                                throw std::runtime_error("empty vector of individuals");
                            }
                        gsl_ran_shuffle(rng.get(), indexes.data(), indexes.size(),
                                        sizeof(std::size_t));
                        rv.emplace(m.first,
                                   move_stack(indexes.size(), std::move(indexes)));
                    }

                return rv;
            }

            template <typename METADATATYPE>
            inline void
            copy_individuals(const std::vector<std::size_t>& buffer,
                             std::int32_t destination,
                             std::vector<METADATATYPE>& metadata)
            // NOTE: the "label" field does not change, in case
            // someone is tracking parents during a simulation.
            {
                for (auto i : buffer)
                    {
                        metadata.push_back(metadata[i]);
                        metadata.back().deme = destination;
                    }
            }

            template <typename METADATATYPE>
            inline void
            apply_copies(const GSLrng_t& rng, const MassMigration& mm,
                         const deme_map_t& deme_map, uint32_t t,
                         std::vector<std::size_t>& buffer,
                         std::vector<METADATATYPE>& metadata)
            {
                auto deme_itr = deme_map.find(mm.source);
                if (deme_itr == end(deme_map))
                    {
                        std::ostringstream o;
                        o << "copies from empty deme " << mm.source << " at time " << t
                          << " attempted";
                        throw DemographyError(o.str());
                    }
                if (mm.fraction < 1.)
                    {
                        std::size_t destination_size = std::round(
                            static_cast<double>(deme_itr->second.size()) * mm.fraction);
                        buffer.resize(destination_size);
                        // Cannot choose an individual 2x to copy to
                        // the same destination.
                        int rv = gsl_ran_choose(
                            rng.get(), buffer.data(), destination_size,
                            const_cast<std::size_t*>(deme_itr->second.data()),
                            deme_itr->second.size(), sizeof(std::size_t));
                        if (rv != GSL_SUCCESS)
                            {
                                throw std::runtime_error(
                                    "MassMigration error: gsl_ran_choose "
                                    "returned "
                                    "failure");
                            }
                        copy_individuals(buffer, mm.destination, metadata);
                    }
                else // entire deme is copied.
                    {
                        copy_individuals(deme_itr->second, mm.destination, metadata);
                    }
            }

            inline void
            move_from_stack(std::int32_t destination, std::size_t n,
                            std::vector<std::int32_t>& moves, move_stack& ms)
            {
                for (std::size_t i = 0; i < n && !ms.ms.empty(); ++i)
                    {
                        auto val = ms.ms.top();
                        if (moves[val] != -1)
                            {
                                throw std::runtime_error(
                                    "MassMigration error: individual has "
                                    "already "
                                    "moved");
                            }
                        moves[val] = destination;
                        ms.ms.pop();
                    }
            }

            inline void
            apply_moves(const MassMigration& mm, std::uint32_t t,
                        std::vector<std::int32_t>& moves, move_map_t& move_source)
            {
                auto deme_itr = move_source.find(mm.source);
                if (deme_itr == end(move_source))
                    {
                        std::ostringstream o;
                        o << "moves from empty deme " << mm.source << " at time " << t
                          << " attempted";
                        throw DemographyError(o.str());
                    }
                if (mm.fraction < 1.0)
                    {
                        std::size_t destination_size
                            = std::round(static_cast<double>(deme_itr->second.initial_N)
                                         * mm.fraction);
                        move_from_stack(mm.destination, destination_size, moves,
                                        deme_itr->second);
                    }
                else // entire deme moves
                    {
                        move_from_stack(mm.destination, deme_itr->second.initial_N,
                                        moves, deme_itr->second);
                    }
            }

            inline void
            update_changed_and_reset(
                const MassMigration& mm,
                std::unordered_map<std::int32_t, bool>& changed_and_reset)
            {
                if (mm.resets_growth_rate == true)
                    {
                        changed_and_reset[mm.source] = true;
                        changed_and_reset[mm.destination] = true;
                    }
                else
                    {
                        auto itr = changed_and_reset.find(mm.source);
                        if (itr == end(changed_and_reset))
                            {
                                changed_and_reset[mm.source] = false;
                            }
                        itr = changed_and_reset.find(mm.destination);
                        if (itr == end(changed_and_reset))
                            {
                                changed_and_reset[mm.destination] = false;
                            }
                    }
            }

            inline std::unordered_map<std::int32_t, std::size_t>
            get_deme_sizes(const deme_map_t& deme_map)
            {
                std::unordered_map<std::int32_t, std::size_t> rv;
                for (auto&& i : deme_map)
                    {
                        rv.emplace(i.first, i.second.size());
                    }
                return rv;
            }

            template <typename METADATATYPE>
            inline std::unordered_map<std::int32_t, std::size_t>
            update_metadata_due_to_moves_and_copies(
                const std::size_t initial_N, const deme_map_t& deme_map,
                const std::vector<std::int32_t>& moves,
                std::vector<METADATATYPE>& metadata)
            {
                auto deme_sizes = detail::get_deme_sizes(deme_map);
                // If there are moves, update the metadata and the deme sizes
                for (std::size_t m = 0; m < moves.size(); ++m)
                    {
                        if (moves[m] != -1)
                            {
                                auto itr = deme_sizes.find(metadata[m].deme);
                                if (itr == deme_sizes.end())
                                    {
                                        throw std::runtime_error(
                                            "MassMigration error: deme not in "
                                            "lookup "
                                            "table");
                                    }
                                if (itr->second == 0)
                                    {
                                        throw std::runtime_error(
                                            "MassMigration error: deme size "
                                            "in "
                                            "lookup "
                                            "table is zero during moves");
                                    }
                                itr->second--;
                                metadata[m].deme = moves[m];
                                deme_sizes[moves[m]]++;
                            }
                    }
                for (std::size_t i = initial_N; i < metadata.size(); ++i)
                    {
                        deme_sizes[metadata[i].deme]++;
                    }
                return deme_sizes;
            }

            inline void
            update_growth_parameters(
                std::uint32_t t,
                const std::unordered_map<std::int32_t, std::size_t>& final_deme_sizes,
                const std::unordered_map<std::int32_t, bool>& changed_and_reset,
                growth_rates_vector& growth_rates,
                growth_rates_onset_times_vector& growth_rate_onset_times,
                growth_initial_size_vector& growth_initial_sizes)
            {
                for (auto&& cr : changed_and_reset)
                    {
                        auto fds = final_deme_sizes.find(cr.first);
                        if (fds == final_deme_sizes.end())
                            {
                                throw std::runtime_error(
                                    "MassMigration error: changed_and_reset "
                                    "key "
                                    "not found in final_deme_sizes");
                            }
                        // NOTE: fds->second is size_t, but deme sizes
                        // are uint32_t so we have to check for overflow
                        if (fds->second >= std::numeric_limits<std::uint32_t>::max())
                            {
                                throw std::runtime_error(
                                    "MassMigration error: deme size overflow");
                            }
                        if (cr.second == true)
                            {
                                growth_rates.get()[cr.first] = NOGROWTH;
                            }
                        growth_rate_onset_times.get()[fds->first] = t;
                        growth_initial_sizes.get()[fds->first] = fds->second;
                    }
            }
        } // namespace detail

        template <typename METADATATYPE>
        inline void
        apply_mass_migrations(const GSLrng_t& rng, std::uint32_t t,
                              current_event_state<MassMigration>& mass_migration_events,
                              growth_rates_vector& growth_rates,
                              growth_rates_onset_times_vector& growth_rate_onset_times,
                              growth_initial_size_vector& growth_initial_sizes,
                              std::vector<METADATATYPE>& metadata)
        // Simultaneous application of all mass migration events
        // occurring at time t
        {
            if (mass_migration_events.current() < mass_migration_events.last()
                && mass_migration_events.when() != t)
                {
                    // TODO: maybe this should be an exception,
                    // or is not necessary at all?
                    return;
                }
            std::vector<METADATATYPE> copies;
            const std::size_t initial_N = metadata.size();
            std::vector<std::int32_t> moves;
            auto deme_map = detail::build_deme_map(metadata);
            std::vector<std::size_t> buffer;
            detail::move_map_t move_source;
            bool initialized_moves = false;

            std::unordered_map<std::int32_t, bool> changed_and_reset;
            for (; mass_migration_events.current() < mass_migration_events.last()
                   && mass_migration_events.when() == t;
                 ++mass_migration_events.current())
                {
                    if (mass_migration_events.event().move_individuals == false) //copy
                        {
                            if (initialized_moves == true)
                                {
                                    // NOTE: this may no longer be necessary?
                                    throw std::runtime_error(
                                        "MassMigration error: copies after "
                                        "moves");
                                }
                            detail::apply_copies(rng, mass_migration_events.event(),
                                                 deme_map, t, buffer, metadata);
                            detail::update_changed_and_reset(
                                mass_migration_events.event(), changed_and_reset);
                        }
                    else //move
                        {
                            if (initialized_moves == false)
                                {
                                    moves.resize(initial_N, -1);
                                    move_source
                                        = detail::build_move_sources(rng, deme_map);
                                    initialized_moves = true;
                                }
                            detail::apply_moves(mass_migration_events.event(), t, moves,
                                                move_source);
                            detail::update_changed_and_reset(
                                mass_migration_events.event(), changed_and_reset);
                        }
                }
            auto final_deme_sizes = detail::update_metadata_due_to_moves_and_copies(
                initial_N, deme_map, moves, metadata);
            detail::update_growth_parameters(t, final_deme_sizes, changed_and_reset,
                                             growth_rates, growth_rate_onset_times,
                                             growth_initial_sizes);
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
