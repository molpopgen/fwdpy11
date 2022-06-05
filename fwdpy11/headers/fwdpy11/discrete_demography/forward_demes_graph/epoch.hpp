#pragma once

#include <cstdint>
#include <memory>
#include "demes_model_time.hpp"
#include "size_function.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct Epoch
        {
            demes_model_time end_time;
            std::uint32_t start_size, end_size;
            double cloning_rate, selfing_rate;
            SizeFunction size_function;

            Epoch(demes_model_time end_time, std::uint32_t start_size,
                  std::uint32_t end_size, double cloning_rate, double selfing_rate,
                  SizeFunction size_function)
                : end_time{end_time}, start_size{start_size}, end_size{end_size},
                  cloning_rate{cloning_rate}, selfing_rate{selfing_rate},
                  size_function{std::move(size_function)}
            {
            }
        };
    }
}
