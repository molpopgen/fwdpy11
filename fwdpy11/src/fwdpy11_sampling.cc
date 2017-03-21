#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/matrix.hpp>
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

    py::class_<KTfwd::data_matrix>(m,"DataMatrix")
        .def(py::init<>())
        .def(py::init<std::size_t>())
        .def_readonly("neutral",&KTfwd::data_matrix::neutral)
        .def_readonly("selected",&KTfwd::data_matrix::selected)
        .def_readonly("neutral_positions",&KTfwd::data_matrix::neutral_positions)
        .def_readonly("selected_positions",&KTfwd::data_matrix::selected_positions)
        .def_readonly("neutral_popfreq",&KTfwd::data_matrix::neutral_popfreq)
        .def_readonly("selected_popfreq",&KTfwd::data_matrix::selected_popfreq)
        .def_readonly("nrow",&KTfwd::data_matrix::nrow)
        ;

    return m.ptr();
}
