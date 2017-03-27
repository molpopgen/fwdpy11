#ifndef FWDPY_SAMPLERS_HPP__
#define FWDPY_SAMPLERS_HPP__

#include <functional>
#include <fwdpy11/types.hpp>

namespace fwdpy11
{
    //Applied each generation to record any data of interest.
    using singlepop_temporal_sampler
        = std::function<void(const fwdpy11::singlepop_t&)>;
}

#endif
