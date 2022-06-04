#pragma once

#include <cstdint>
#include "demes_model_time.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct AsymmetricMigration
        {
            std::int32_t source, dest;
            demes_model_time start_time, end_time;
        };
    }
}
