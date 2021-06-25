#pragma once

#include <gsl/gsl_randist.h>
#include <stack>
#include <vector>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <fwdpy11/discrete_demography/DiscreteDemographyState.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace fwdpy11
{
    namespace discrete_demography
    {
        class DiscreteDemographyState::DiscreteDemographyState_impl
        {
          private:
            template <typename T> struct events_with_range
            {
                std::vector<T> events;
                std::pair<std::size_t, std::size_t> event_range;
                template <typename Input>
                events_with_range(Input&& input)
                    : events(std::forward<Input>(input)), event_range{0, events.size()}
                {
                }
            };

            // These are the input parameters to the model +
            // a set of integers to track what the next event is.
            // The event lists and the migration matrix M are de-facto
            // const.
            events_with_range<MassMigration> mass_migrations;
            events_with_range<SetExponentialGrowth> set_growth_rates;
            events_with_range<SetDemeSize> set_deme_sizes;
            events_with_range<SetSelfingRate> set_selfing_rates;
            MigrationMatrix M;
            events_with_range<SetMigrationRates> set_migration_rates;

            // The current state of the model
            std::vector<std::uint32_t> current_deme_sizes;
            std::vector<std::uint32_t> next_deme_sizes;
            std::vector<std::uint32_t> growth_onset_times;
            std::vector<std::uint32_t> growth_initial_sizes;
            std::vector<double> growth_rates;
            std::vector<double> selfing_rates;
            std::uint32_t current_time_in_simulation;
            std::uint32_t next_global_N;

            // Implementation details of demographic events
            // Template implementations are in case of future
            // support for non-diploid populations.
            // Much of these details are copied and slightly modified
            // from standalone functions prior to 0.16.0.

            using deme_map = std::unordered_map<std::int32_t, std::vector<std::size_t>>;
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

            using move_map = std::unordered_map<std::int32_t, move_stack>;

            std::unordered_map<std::int32_t, std::size_t>
            get_deme_sizes(const deme_map& deme_map)
            {
                std::unordered_map<std::int32_t, std::size_t> rv;
                for (auto&& i : deme_map)
                    {
                        rv.emplace(i.first, i.second.size());
                    }
                return rv;
            }

            template <typename METADATATYPE>
            deme_map
            build_deme_map(const std::vector<METADATATYPE>& metadata)
            {
                deme_map rv;
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

            move_map
            build_move_sources(const GSLrng_t& rng, const deme_map& deme_map)
            {
                move_map rv;

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

            // TODO: better name for this
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

            template <typename METADATATYPE>
            inline void
            mass_migration_copy_individuals(const std::vector<std::size_t>& buffer,
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

            void
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
            mass_migration_moves(const MassMigration& mm, std::uint32_t t,
                                 std::vector<std::int32_t>& moves, move_map& move_source)
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

            template <typename METADATATYPE>
            inline void
            mass_migration_copies(const GSLrng_t& rng, const MassMigration& mm,
                                  const deme_map& deme_map, uint32_t t,
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
                        mass_migration_copy_individuals(buffer, mm.destination,
                                                        metadata);
                    }
                else // entire deme is copied.
                    {
                        mass_migration_copy_individuals(deme_itr->second, mm.destination,
                                                        metadata);
                    }
            }

            template <typename METADATATYPE>
            std::unordered_map<std::int32_t, std::size_t>
            update_metadata_due_to_mass_migrations(
                const std::size_t initial_N, const deme_map& deme_map,
                const std::vector<std::int32_t>& moves,
                std::vector<METADATATYPE>& metadata)
            {
                auto deme_sizes = get_deme_sizes(deme_map);
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

            void
            update_growth_parameters_after_mass_migrations(
                std::uint32_t t,
                const std::unordered_map<std::int32_t, std::size_t>& final_deme_sizes,
                const std::unordered_map<std::int32_t, bool>& changed_and_reset)
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
                                growth_rates[cr.first] = NOGROWTH;
                            }
                        growth_onset_times[fds->first] = t;
                        growth_initial_sizes[fds->first] = fds->second;
                    }
            }

            template <typename METADATATYPE>
            void
            apply_mass_migrations(const GSLrng_t& rng,
                                  const std::uint32_t simulation_time,
                                  std::vector<METADATATYPE>& individual_metadata)
            {
                if (mass_migrations.event_range.first
                        < mass_migrations.event_range.second
                    && mass_migrations.events[mass_migrations.event_range.first].when
                           != simulation_time)
                    {
                        return;
                    }
                for (; mass_migrations.event_range.first
                           < mass_migrations.event_range.second
                       && mass_migrations.events[mass_migrations.event_range.first].when
                              == simulation_time;
                     ++mass_migrations.event_range.first)
                    {
                        const auto initial_N = individual_metadata.size();
                        if (initial_N == 0)
                            {
                                throw std::runtime_error("metadata are empty");
                            }
                        bool initialized_moves{false};
                        std::vector<std::int32_t> moves; // TODO: is this type okay?
                        std::vector<std::size_t> buffer;

                        auto individual_to_deme = build_deme_map(individual_metadata);
                        move_map move_source;
                        std::unordered_map<std::int32_t, bool> changed_and_reset;

                        for (std::size_t i = mass_migrations.event_range.first;
                             i < mass_migrations.event_range.second
                             && mass_migrations.events[i].when == simulation_time;
                             ++i)
                            {
                                if (mass_migrations.events[i].move_individuals
                                    == false) // copy
                                    {
                                        if (initialized_moves == true)
                                            {
                                                // NOTE: this may no longer be necessary?
                                                throw std::runtime_error(
                                                    "MassMigration error: copies after "
                                                    "moves");
                                            }
                                        mass_migration_copies(
                                            rng, mass_migrations.events[i],
                                            individual_to_deme, simulation_time, buffer,
                                            individual_metadata);
                                    }
                                else // move event
                                    {
                                        if (initialized_moves == false)
                                            {
                                                moves.resize(initial_N, -1);
                                                move_source = build_move_sources(
                                                    rng, individual_to_deme);
                                                initialized_moves = true;
                                            }
                                        mass_migration_moves(mass_migrations.events[i],
                                                             simulation_time, moves,
                                                             move_source);
                                    }
                                update_changed_and_reset(mass_migrations.events[i],
                                                         changed_and_reset);
                            }
                        auto final_deme_sizes = update_metadata_due_to_mass_migrations(
                            initial_N, individual_to_deme, moves, individual_metadata);
                        update_growth_parameters_after_mass_migrations(
                            simulation_time, final_deme_sizes, changed_and_reset);
                    }
            }

          public:
            DiscreteDemographyState_impl(
                std::vector<MassMigration> mass_migrations,
                std::vector<SetExponentialGrowth> set_growth_rates,
                std::vector<SetDemeSize> set_deme_sizes,
                std::vector<SetSelfingRate> set_selfing_rates, MigrationMatrix M,
                std::vector<SetMigrationRates> set_migration_rates);

            // This function must be run early, when the simulation
            // code is entered, but before we start iterating over generations.
            // This function needs to take some/all responsibility for
            // determining the max number of demes in the simulation.
            void initialize(const DiploidPopulation& pop,
                            const std::uint32_t pop_generation);

            // Applies mass migration events and deme size changes.
            // Will affect growh rates, too.
            // NOTE: replaces mass_migrations_and_current_sizes.
            template <typename METADATATYPE>
            void
            early(const GSLrng_t& rng, const std::uint32_t generation,
                  std::vector<METADATATYPE>& individual_metadata)
            {
                apply_mass_migrations(rng, generation, individual_metadata);
            }

            // Updates fitness lookups, migration lookups,
            // and performs runtime validations.
            // NOTE: replaces finalize_demographic_state.
            template <typename METADATATYPE>
            void late(const GSLrng_t& rng, const std::uint32_t generation,
                      std::vector<METADATATYPE>& individual_metadata);

            // Assist the DiscreteDemographyState copy constructor
            std::unique_ptr<DiscreteDemographyState_impl> clone() const;
        };

        extern template void
        DiscreteDemographyState::DiscreteDemographyState_impl::early(
            const GSLrng_t& rng, const std::uint32_t generation,
            std::vector<DiploidMetadata>& metadata);

        extern template void DiscreteDemographyState::DiscreteDemographyState_impl::late(
            const GSLrng_t& rng, const std::uint32_t generation,
            std::vector<DiploidMetadata>& metadata);
    }
}
