#include "DiscreteDemographyState_impl.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        DiscreteDemographyState::DiscreteDemographyState_impl::
            DiscreteDemographyState_impl(
                std::vector<MassMigration> mass_migrations,
                std::vector<SetExponentialGrowth> set_growth_rates,
                std::vector<SetDemeSize> set_deme_sizes,
                std::vector<SetSelfingRate> set_selfing_rates, MigrationMatrix M,
                std::vector<SetMigrationRates> set_migration_rates)
            : mass_migrations{std::move(mass_migrations)},
              set_growth_rates{std::move(set_growth_rates)}, set_deme_sizes{std::move(
                                                                 set_deme_sizes)},
              set_selfing_rates{std::move(set_selfing_rates)}, M{std::move(M)},
              set_migration_rates{std::move(set_migration_rates)}, current_deme_sizes{},
              next_deme_sizes{}, growth_onset_times{}, growth_initial_sizes{},
              growth_rates{}, selfing_rates{}, current_time_in_simulation{0}
        {
        }

        std::unique_ptr<DiscreteDemographyState::DiscreteDemographyState_impl>
        DiscreteDemographyState::DiscreteDemographyState_impl::clone() const
        {
            return std::make_unique<
                DiscreteDemographyState::DiscreteDemographyState_impl>(*this);
        }

        void
        DiscreteDemographyState::DiscreteDemographyState_impl::initialize(
            const std::vector<std::uint32_t>& deme_sizes,
            const std::uint32_t pop_generation)
        {
        }
    }
}
