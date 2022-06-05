#pragma once

#include <functional>
#include <stdexcept>
#include "demes_model_time.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        // Should this be a virtual base class
        // or simply store a std::function?
        struct SizeFunction
        {
            using size_function = std::function<std::uint32_t(
                std::uint32_t /*epoch_start_size*/, std::uint32_t /*epoch_end_size*/,
                demes_model_time /*epoch_start_time*/,
                demes_model_time /*epoch_end_time*/, demes_model_time /*current_time*/)>;
            size_function f;

            SizeFunction(size_function f) : f{f}
            {
            }
        };

        struct ConstantSizeFunction
        {
            inline std::uint32_t
            operator()(std::uint32_t epoch_start_size, std::uint32_t epoch_end_size,
                       demes_model_time /*epoch_start_time*/,
                       demes_model_time /*epoch_end_time*/,
                       demes_model_time /*current_time*/)
            {
                if (epoch_start_size != epoch_end_size)
                    {
                        throw std::invalid_argument(
                            "start_size != end_size incompatible with "
                            "ConstantSizeFunction");
                    }
                return epoch_start_size;
            }
        };
    }
}
