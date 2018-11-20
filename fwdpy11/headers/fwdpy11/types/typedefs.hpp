#ifndef FWDPY11_TYPES_TYPEDEFS_HPP__
#define FWDPY11_TYPES_TYPEDEFS_HPP__

#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include "Diploid.hpp"

namespace fwdpy11
{
    //! Typedef for mutation container
    using mcont_t = std::vector<fwdpp::popgenmut>;
    //! Typedef for gamete type
    using gamete_t = fwdpp::gamete;
    //! Typedef for gamete container
    using gcont_t = std::vector<gamete_t>;
}
#endif
