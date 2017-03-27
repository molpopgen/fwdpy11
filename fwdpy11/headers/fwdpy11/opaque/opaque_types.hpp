#ifndef FWDPY11_OPAQUE_TYPES_HPP__
#define FWDPY11_OPAQUE_TYPES_HPP__

#include <pybind11/stl.h>
#include <cstdint>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace fwdpy11
{
    //! Typedef for mutation container
    using mcont_t = std::vector<KTfwd::popgenmut>;
    //! Typedef for gamete type
    using gamete_t = KTfwd::gamete;
    //! Typedef for gamete container
    using gcont_t = std::vector<gamete_t>;
}

PYBIND11_MAKE_OPAQUE(fwdpy11::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::mcont_t);

#endif
