#pragma once

#include <functional>
#include <stdexcept>
#include "demes_model_time.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct SizeFunction
        {
            using size_function = std::function<std::uint32_t(
                std::uint32_t /*epoch_start_size*/, std::uint32_t /*epoch_end_size*/,
                demes_model_time /*epoch_start_time*/,
                demes_model_time /*epoch_end_time*/, demes_model_time /*current_time*/)>;

            // Throw std::invalid_argument if start/end sizes
            // are not compatible
            using validate_function = std::function<void(
                std::uint32_t /*epoch_start_size*/, std::uint32_t /*epoch_end_size*/)>;
            size_function size_at;
            validate_function validate;

            SizeFunction(size_function size_at, validate_function validate)
                : size_at{std::move(size_at)}, validate{std::move(validate)}
            {
            }

            inline static SizeFunction
            constant()
            {
                return SizeFunction{
                    [](std::uint32_t epoch_start_size, std::uint32_t /*epoch_end_size*/,
                       demes_model_time /*epoch_start_time*/,
                       demes_model_time /*epoch_end_time*/,
                       demes_model_time /*current_time*/) { return epoch_start_size; },
                    [](std::uint32_t epoch_start_size, std::uint32_t epoch_end_size) {
                        if (epoch_start_size != epoch_end_size)
                            {
                                throw std::invalid_argument(
                                    "start_size != end_size incompatible with constant "
                                    "size function");
                            }
                    }};
            }
        };
    }
}
