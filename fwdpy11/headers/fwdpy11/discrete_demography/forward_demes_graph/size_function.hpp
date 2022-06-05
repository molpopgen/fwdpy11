#pragma once

#include <functional>
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
                demes_model_time /*epoch_start_time*/,
                demes_model_time /*epoch_end_time*/, demes_model_time /*current_time*/)>;
            size_function f;

            SizeFunction(size_function f) : f{f}
            {
            }
        };

    }
}
