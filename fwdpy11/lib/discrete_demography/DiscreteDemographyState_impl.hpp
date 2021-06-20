#pragma once

#include <gsl/gsl_randist.h>
#include <stack>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <fwdpy11/discrete_demography/DiscreteDemographyState.hpp>

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
                const auto initial_N = individual_metadata.size();
                if (initial_N == 0)
                    {
                        throw std::runtime_error("metadata are empty");
                    }
                bool initialized_moves{false};
                std::vector<METADATATYPE> copies;
                std::vector<std::int32_t> moves; // TODO: is this type okay?
                std::vector<std::size_t> buffer;

                auto input_deme_map = build_deme_map(individual_metadata);
                move_map move_source;

                for (std::size_t i = mass_migrations.event_range.first;
                     i < mass_migrations.event_range.second
                     && mass_migrations.events[i].when == simulation_time;
                     ++i)
                    {
                        if (mass_migrations.events[i].move_individuals == false) // copy
                            {
                                if (initialized_moves == true)
                                    {
                                        // NOTE: this may no longer be necessary?
                                        throw std::runtime_error(
                                            "MassMigration error: copies after "
                                            "moves");
                                    }
                            }
                        else // move event
                            {
                            }
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
            void initialize(const std::vector<std::uint32_t>& deme_sizes,
                            const std::uint32_t pop_generation);

            // Applies mass migration events and deme size changes.
            // Will affect growh rates, too.
            // NOTE: replaces mass_migrations_and_current_sizes.
            void early(const GSLrng_t& rng, const std::uint32_t generation,
                       std::vector<DiploidMetadata>& metadata);

            // Updates fitness lookups, migration lookups,
            // and performs runtime validations.
            // NOTE: replaces finalize_demographic_state.
            void late(const GSLrng_t& rng, const std::uint32_t generation,
                      std::vector<DiploidMetadata>& metadata);

            // Assist the DiscreteDemographyState copy constructor
            std::unique_ptr<DiscreteDemographyState_impl> clone() const;
        };
    }
}
