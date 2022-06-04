#pragma once

#include <cstdint>

namespace fwdpy11
{
    namespace discrete_demography
    {
        // fwdpy11 time is uint32_t, but that is
        // BAD if you have to do stuff like subtract!
        // So, we use int64_t here and will return casted
        // values as needed (after checks).
        using demes_model_time = std::int64_t;
    }
}
