#include <utility>
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

            events_with_range<MassMigration> mass_migrations;
            events_with_range<SetExponentialGrowth> set_growth_rates;
            events_with_range<SetDemeSize> set_deme_sizes;
            events_with_range<SetSelfingRate> set_selfing_rates;
            MigrationMatrix M;
            events_with_range<SetMigrationRates> set_migration_rates;

          public:
            DiscreteDemographyState_impl(
                std::vector<MassMigration> mass_migrations,
                std::vector<SetExponentialGrowth> set_growth_rates,
                std::vector<SetDemeSize> set_deme_sizes,
                std::vector<SetSelfingRate> set_selfing_rates, MigrationMatrix M,
                std::vector<SetMigrationRates> set_migration_rates)
                : mass_migrations{std::move(mass_migrations)},
                  set_growth_rates{std::move(set_growth_rates)},
                  set_deme_sizes{std::move(set_deme_sizes)},
                  set_selfing_rates{std::move(set_selfing_rates)}, M{std::move(M)},
                  set_migration_rates{std::move(set_migration_rates)}
            {
            }

            std::unique_ptr<DiscreteDemographyState_impl>
            clone() const
            {
                return std::make_unique<DiscreteDemographyState_impl>(*this);
            }
        };

        DiscreteDemographyState::DiscreteDemographyState(
            std::vector<MassMigration> mass_migrations,
            std::vector<SetExponentialGrowth> set_growth_rates,
            std::vector<SetDemeSize> set_deme_sizes,
            std::vector<SetSelfingRate> set_selfing_rates, MigrationMatrix M,
            std::vector<SetMigrationRates> set_migration_rates)
            : pimpl{new DiscreteDemographyState_impl{
                std::move(mass_migrations), std::move(set_growth_rates),
                std::move(set_deme_sizes), std::move(set_selfing_rates), std::move(M),
                std::move(set_migration_rates)}}
        {
        }
    }
}
