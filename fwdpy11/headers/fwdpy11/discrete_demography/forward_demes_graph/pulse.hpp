#pragma once

#include <vector>
#include <cstdint>
#include "demes_model_time.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct Pulse
        {
            std::vector<std::int32_t> sources;
            std::vector<double> proportions;
            std::int32_t dest;
            demes_model_time time;
        };
    }
}
