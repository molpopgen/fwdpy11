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
            virtual std::uint32_t size_at(std::uint32_t /*epoch_start_size*/,
                                          std::uint32_t /*epoch_end_size*/,
                                          demes_model_time /*epoch_start_time*/,
                                          demes_model_time /*epoch_end_time*/,
                                          demes_model_time /*current_time*/) const = 0;

            // Throw std::invalid_argument if start/end sizes
            // are not compatible
            virtual void validate(std::uint32_t /*epoch_start_size*/,
                                  std::uint32_t /*epoch_end_size*/) const = 0;
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
