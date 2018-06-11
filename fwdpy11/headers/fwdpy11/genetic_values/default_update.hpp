#ifndef FWDPY11_GENETIC_VALUE_DEFAULT_UPDATE_HPP__
#define FWDPY11_GENETIC_VALUE_DEFAULT_UPDATE_HPP__

#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

#define DEFAULT_SLOCUSPOP_UPDATE()                                            \
    inline void update(const fwdpy11::SlocusPop& /*pop*/) {}

#define DEFAULT_MLOCUSPOP_UPDATE()                                            \
    inline void update(const fwdpy11::MlocusPop& /*pop*/) {}

#endif
