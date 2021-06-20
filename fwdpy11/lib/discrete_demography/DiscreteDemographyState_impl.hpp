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
            // a set of integers to track what the next event is
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

          public:
            DiscreteDemographyState_impl(
                std::vector<MassMigration> mass_migrations,
                std::vector<SetExponentialGrowth> set_growth_rates,
                std::vector<SetDemeSize> set_deme_sizes,
                std::vector<SetSelfingRate> set_selfing_rates, MigrationMatrix M,
                std::vector<SetMigrationRates> set_migration_rates);

            std::unique_ptr<DiscreteDemographyState_impl> clone() const;
        };
    }
}
