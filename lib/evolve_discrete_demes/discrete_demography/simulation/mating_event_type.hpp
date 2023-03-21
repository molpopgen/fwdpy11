#ifndef FWDPY11_DISCRETE_DEMOGRAPY_MATING_EVENT_TYPE_HPP
#define FWDPY11_DISCRETE_DEMOGRAPY_MATING_EVENT_TYPE_HPP

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        enum class mating_event_type
        {
            outcrossing,
            selfing,
            cloning // NOTE: not supported
        };
    }
}

#endif
