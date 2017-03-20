#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include "types.hpp"

namespace py = pybind11;

using sample_t = std::vector<std::pair<double,std::string>>;
using sep_sample_t = std::pair<sample_t,sample_t>;

sep_sample_t sample_singlepop(const fwdpy::GSLrng_t & rng, const fwdpy::singlepop_t & pop,
        const unsigned nsam, const bool removeFixed)
{
    return KTfwd::sample_separate(rng.get(),pop,nsam,removeFixed);
}

PYBIND11_PLUGIN(fwdpy11_sampling) {
    py::module m("fwdpy11_sampling", "Taking samples from populations");

    m.def("sample_singlepop",&sample_singlepop);

    return m.ptr();
}
