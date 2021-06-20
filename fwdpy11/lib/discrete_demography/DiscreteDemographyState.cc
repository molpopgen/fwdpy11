#include <utility>
#include <fwdpy11/discrete_demography/DiscreteDemographyState.hpp>
#include "DiscreteDemographyState_impl.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
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

        DiscreteDemographyState::DiscreteDemographyState(
            const DiscreteDemographyState& other)
            : pimpl{other.pimpl->clone()}
        {
        }

        DiscreteDemographyState::DiscreteDemographyState(DiscreteDemographyState&& other)
            : pimpl{std::exchange(other.pimpl, nullptr)}
        {
        }

        void
        DiscreteDemographyState::early(const GSLrng_t& rng,
                                       const std::uint32_t generation,
                                       std::vector<DiploidMetadata>& metadata)
        {
            pimpl->early(rng, generation, metadata);
        }

        void
        DiscreteDemographyState::late(const GSLrng_t& rng,
                                      const std::uint32_t generation,
                                      std::vector<DiploidMetadata>& metadata)
        {
            pimpl->late(rng, generation, metadata);
        }

    }
}
