#ifndef FWDPY11_GENETIC_VALUES_SLOCUSMULT_HPP__
#define FWDPY11_GENETIC_VALUES_SLOCUSMULT_HPP__

#include <fwdpp/fitness_models.hpp>
#include "fwdpp_wrappers/fwdpp_slocus_gvalue.hpp"

namespace fwdpy11
{
    using SlocusMult = fwdpp_slocus_gvalue<fwdpp::multiplicative_diploid>;
} // namespace fwdpy11
#endif
