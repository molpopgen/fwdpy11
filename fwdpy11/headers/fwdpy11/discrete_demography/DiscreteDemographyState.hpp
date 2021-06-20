#ifndef FWDPY11_DISCRETE_DEMOGRAPHY_STATE_HPP
#define FWDPY11_DISCRETE_DEMOGRAPHY_STATE_HPP

#include <memory>
#include <vector>
#include <utility>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpy11/rng.hpp>
#include "MassMigration.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetDemeSize.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "MigrationMatrix.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class DiscreteDemographyState
        {
          private:
            class DiscreteDemographyState_impl;
            std::unique_ptr<DiscreteDemographyState_impl> pimpl;

          public:
            DiscreteDemographyState(std::vector<MassMigration> mass_migrations,
                                    std::vector<SetExponentialGrowth> set_growth_rates,
                                    std::vector<SetDemeSize> size_changes,
                                    std::vector<SetSelfingRate> set_selfing_rates,
                                    MigrationMatrix M,
                                    std::vector<SetMigrationRates> set_migration_rates);

            DiscreteDemographyState(const DiscreteDemographyState&);
            DiscreteDemographyState(DiscreteDemographyState&&);

            // Applies mass migration events and deme size changes.
            // Will affect growh rates, too.
            void early(const GSLrng_t& rng, const std::uint32_t generation,
                       std::vector<DiploidMetadata>& metadata);

            // Updates fitness lookups, migration lookups,
            // and performs runtime validations.
            void late(const GSLrng_t& rng, const std::uint32_t generation,
                      std::vector<DiploidMetadata>& metadata);
        };

    }
}

#endif

