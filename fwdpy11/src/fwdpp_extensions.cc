#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

namespace py = pybind11;

using dfe_callback_type = std::function<double(const gsl_rng *)>;

template <typename S> struct make_sh_model_fixed_dom
{
    template <typename... args>
    KTfwd::extensions::shmodel
    operator()(const double h, args... A) const
    {
        return KTfwd::extensions::shmodel(
            std::bind(S(A...), std::placeholders::_1),
            std::bind(KTfwd::extensions::constant(h), std::placeholders::_1));
    }
};

#define RETURN_DFE_FIXEDH(DIST, S, H)                                         \
    return make_sh_model_fixed_dom<KTfwd::extensions::DIST>()(S, H);
#define RETURN_DFE2_FIXEDH(DIST, A, B, H)                                     \
    return make_sh_model_fixed_dom<KTfwd::extensions::DIST>()(A, B, H);

PYBIND11_PLUGIN(fwdpp_extensions)
{
    py::module m("fwdpp_extensions", "Expose fwdpp's extensions library.");

    py::class_<KTfwd::extensions::shmodel>(m, "DFEFixedDominance")
        .def(py::init<>())
        .def(py::init<dfe_callback_type, dfe_callback_type>());

    m.def("makeConstantSH", ([](const double s, const double h) {
              RETURN_DFE_FIXEDH(constant, s, h);
          }));

    m.def("makeExpSH", ([](const double mean, const double h) {
              RETURN_DFE_FIXEDH(exponential, mean, h);
          }));

    m.def("makeGaussianSH", ([](const double sd, const double h) {
              RETURN_DFE_FIXEDH(gaussian, sd, h);
          }));

    m.def("makeUniformSH",
          ([](const double lo, const double hi, const double h) {
              RETURN_DFE2_FIXEDH(uniform, lo, hi, h);
          }));

    m.def("makeGammaSH",
          ([](const double mean, const double shape, const double h) {
              RETURN_DFE2_FIXEDH(gamma, mean, shape, h);
          }));
	
    py::class_<KTfwd::extensions::discrete_mut_model>(m, "MutationRegions")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<KTfwd::extensions::shmodel>>())
        .def_readonly("nbeg", &KTfwd::extensions::discrete_mut_model::nbeg)
        .def_readonly("nend", &KTfwd::extensions::discrete_mut_model::nend)
        .def_readonly("sbeg", &KTfwd::extensions::discrete_mut_model::sbeg)
        .def_readonly("send", &KTfwd::extensions::discrete_mut_model::send)
        .def_readonly("shmodels",
                      &KTfwd::extensions::discrete_mut_model::shmodels);

    py::class_<KTfwd::extensions::discrete_rec_model>(m,
                                                      "RecombinationRegions")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>>());

    return m.ptr();
}
