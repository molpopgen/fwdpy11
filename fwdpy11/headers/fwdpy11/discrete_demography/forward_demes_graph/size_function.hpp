#pragma once

#include "demes_model_time.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        // Should this be a virtual base class
        // or simply store a std::function?
        struct SizeFunction
        {
            virtual double operator()(demes_model_time epoch_start_time,
                                      demes_model_time epoch_end_time,
                                      demes_model_time current_time) const;
        };

    }
}
