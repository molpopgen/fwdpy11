#pragma once

#include <cstdint>
#include <memory>
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
            std::vector<Epoch> epochs;

            Deme(std::int32_t id, demes_model_time start_time)
                : id{id}, start_time{start_time}, ancestors{}, proportions{}, epochs{}
            {
            }

            void
            add_epoch(demes_model_time end_time, std::uint32_t start_size,
                      std::uint32_t end_size, double cloning_rate, Selfing selfing,
                      std::unique_ptr<SizeFunction> size_function)
            {
                epochs.push_back(Epoch{end_time, start_size, end_size, cloning_rate,
                                       std::move(selfing), std::move(size_function)});
            }
        };
    }
}
