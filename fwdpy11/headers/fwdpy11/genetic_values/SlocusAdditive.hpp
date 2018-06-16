#ifndef FWDPY11_GENETIC_VALUES_SLOCUSADDITIVE_HPP__
#define FWDPY11_GENETIC_VALUES_SLOCUSADDITIVE_HPP__

#include <fwdpp/fitness_models.hpp>
#include "fwdpp_wrappers/fwdpp_slocus_gvalue.hpp"

namespace fwdpy11
{
    using SlocusAdditive = fwdpp_slocus_gvalue<fwdpp::additive_diploid>;
} // namespace fwdpy11

#endif
