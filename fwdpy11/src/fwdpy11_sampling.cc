#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpy11/types.hpp>

namespace py = pybind11;

template <typename poptype>
typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                     KTfwd::sugar::SINGLEPOP_TAG>::value,
                        KTfwd::sep_sample_t>::type
sample_separate_wrapper(const fwdpy11::GSLrng_t &rng, const poptype &pop,
                        const unsigned nsam, const bool removeFixed)
{
    return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed);
}

template <typename poptype>
typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                     KTfwd::sugar::MULTILOCPOP_TAG>::value,
                        std::vector<KTfwd::sep_sample_t>>::type
sample_separate_wrapper(
    const fwdpy11::GSLrng_t &rng, const poptype &pop, const unsigned nsam,
    const bool removeFixed,
    const std::vector<std::pair<double, double>> &locus_boundaries)
{
    return KTfwd::sample_separate(rng.get(), pop, nsam, removeFixed,
                                  locus_boundaries);
}

PYBIND11_PLUGIN(sampling)
{
    py::module m("sampling", "Taking samples from populations");

    m.def("sample_separate", &sample_separate_wrapper<fwdpy11::singlepop_t>);
    m.def("sample_separate", &sample_separate_wrapper<fwdpy11::multilocus_t>);

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
        .def_readonly("selected_popfreq",
                      &KTfwd::data_matrix::selected_popfreq)
        .def_readonly("nrow", &KTfwd::data_matrix::nrow)
        .def("__getstate__",
             [](const KTfwd::data_matrix &d) {
                 return py::make_tuple(d.nrow, d.neutral, d.selected,
                                       d.neutral_positions,
                                       d.selected_positions, d.neutral_popfreq,
                                       d.selected_popfreq);
             })
        .def("__setstate__", [](KTfwd::data_matrix &d, py::tuple p) {
            new (&d) KTfwd::data_matrix(p[0].cast<std::size_t>());
            d.nrow = p[0].cast<std::size_t>();
            d.neutral = p[1].cast<std::vector<char>>();
            d.selected = p[1].cast<std::vector<char>>();
            d.neutral_positions = p[1].cast<std::vector<double>>();
            d.selected_positions = p[1].cast<std::vector<double>>();
            d.neutral_popfreq = p[1].cast<std::vector<double>>();
            d.selected_popfreq = p[1].cast<std::vector<double>>();
        });

    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy11::singlepop_t>);
    m.def("mutation_keys", &KTfwd::mutation_keys<fwdpy11::multilocus_t>);
    m.def("genotype_matrix", &KTfwd::genotype_matrix<fwdpy11::singlepop_t>);
    m.def("genotype_matrix", &KTfwd::genotype_matrix<fwdpy11::multilocus_t>);
    return m.ptr();
}
