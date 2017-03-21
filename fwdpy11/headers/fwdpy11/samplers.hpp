#ifndef FWDPY_SAMPLERS_HPP__
#define FWDPY_SAMPLERS_HPP__

#include <functional>
#include <fwdpy11/types.hpp>

namespace fwdpy11
{
	//! Takes a const reference to a population
	// and the generation as arguments.
    using singlepop_temporal_sampler
        = std::function<void(const fwdpy11::singlepop_t&, const unsigned)>;
}

#endif
