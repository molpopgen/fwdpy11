#pragma once

#include <fwdpy11/rng.hpp>
#include <gsl/gsl_randist.h>

namespace fwdpy11_core
{
    namespace internal
    {
        inline auto
        gsl_ran_flat(const fwdpy11::GSLrng_t& rng, double lo, double hi)
        {
            auto rv = gsl_ran_flat(rng.get(), lo, hi);
            while (rv == hi)
                {
                    rv = gsl_ran_flat(rng.get(), lo, hi);
                }
            return rv;
        }
    }
}
