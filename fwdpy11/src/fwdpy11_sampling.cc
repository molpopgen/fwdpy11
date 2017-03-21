#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include "types.hpp"

namespace py = pybind11;

template <typename poptype>
typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                     KTfwd::sugar::SINGLEPOP_TAG>::value,
                        KTfwd::sep_sample_t>::type
sample_separate_wrapper(const fwdpy::GSLrng_t &rng, const poptype &pop,
                        const unsigned nsam, const bool removeFixed) {
    return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed);
}

template <typename poptype>
typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                     KTfwd::sugar::MULTILOCPOP_TAG>::value,
                        std::vector<KTfwd::sep_sample_t>>::type
sample_separate_wrapper(
    const fwdpy::GSLrng_t &rng, const poptype &pop, const unsigned nsam,
    const bool removeFixed,
    const std::vector<std::pair<double, double>> &locus_boundaries) {
    return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed,
                                  locus_boundaries);
}

PYBIND11_PLUGIN(fwdpy11_sampling) {
    py::module m("fwdpy11_sampling", "Taking samples from populations");

    m.def("sample_separate", &sample_separate_wrapper<fwdpy::singlepop_t>);
    //m.def("sample_separate", &sample_separate_wrapper<fwdpy::multilocus_t>);

    py::class_<KTfwd::data_matrix>(m, "DataMatrix")
        .def(py::init<>())
        .def(py::init<std::size_t>())
        .def_readonly("neutral", &KTfwd::data_matrix::neutral)
        .def_readonly("selected", &KTfwd::data_matrix::selected)
        .def_readonly("neutral_positions",
                      &KTfwd::data_matrix::neutral_positions)
        .def_readonly("selected_positions",
                      &KTfwd::data_matrix::selected_positions)
        .def_readonly("neutral_popfreq", &KTfwd::data_matrix::neutral_popfreq)
        .def_readonly("selected_popfreq", &KTfwd::data_matrix::selected_popfreq)
        .def_readonly("nrow", &KTfwd::data_matrix::nrow);

    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy::singlepop_t>);
    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy::multilocus_t>);

    return m.ptr();
}
