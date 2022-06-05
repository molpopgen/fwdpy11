#pragma once

#include <cstdint>
#include <vector>
#include "demes_model_time.hpp"
#include "epoch.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct Deme
        {
            std::int32_t id;
            demes_model_time start_time;
            std::vector<std::int32_t> ancestors;
            std::vector<double> proportions;
            //std::vector<Epoch> epochs;

            Deme(std::int32_t id, demes_model_time start_time)
                : id{id}, start_time{start_time}, ancestors{}, proportions{}//, epochs{}
            {
            }
        };
    }
}
