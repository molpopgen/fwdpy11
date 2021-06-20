#pragma once

#include <vector>
#include <cstdint>
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

            template <typename METADATATYPE>
            void
            apply_mass_migrations(const GSLrng_t& rng,
                                  const std::uint32_t simulation_time,
                                  std::vector<METADATATYPE>& individual_metadata)
            {
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
