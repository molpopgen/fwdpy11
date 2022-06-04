#pragma once

#include <memory>
#include "demes_model_time.hpp"
#include "size_function.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct Epoch
        {
            demes_model_time end_time, start_size, end_size;
            double cloning_rate, selfing_rate;

            // Never accept nullptr: "constant" should
            // also be a callback.
            std::unique_ptr<SizeFunction> size_function;
        };
    }
}
