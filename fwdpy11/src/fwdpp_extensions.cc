#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

namespace py = pybind11;

using dfe_callback_type = std::function<double(const gsl_rng *)>;

PYBIND11_PLUGIN(fwdpp_extensions)
{
    py::module m("fwdpp_extensions", "Expose fwdpp's extensions library.");

    py::class_<KTfwd::extensions::shmodel>(m, "DFEFixedDominance")
        .def(py::init<>())
        .def(py::init<dfe_callback_type, dfe_callback_type>());

#define DFE(A, B, C)                                                          \
    py::class_<KTfwd::extensions::A>(m, B).def(C).def_property_readonly(      \
        "callback", [](const KTfwd::extensions::A &c) -> dfe_callback_type {  \
            return std::bind(c, std::placeholders::_1);                       \
        });

    DFE(constant, "ConstantSH", py::init<double>());
    DFE(exponential, "ExpSH", py::init<double>());
    DFE(gaussian, "GaussianSH", py::init<double>());
    DFE(uniform, "UniformSH", (py::init<double, double>()));
    DFE(gamma, "GammaSH", (py::init<double, double>()));
    DFE(beta, "BetaSH", (py::init<double, double>()));

    // py::class_<KTfwd::extensions::constant>(m, "ConstantSH")
    //    .def(py::init<double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::constant &c) -> dfe_callback_type {
    //            return std::bind(c, std::placeholders::_1);
    //        });

    // py::class_<KTfwd::extensions::exponential>(m, "ExpSH")
    //    .def(py::init<double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::exponential &c) -> dfe_callback_type
    //        {
    //            return std::bind(c, std::placeholders::_1);
    //        });

    // py::class_<KTfwd::extensions::uniform>(m, "UniformSH")
    //    .def(py::init<double,double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::uniform &c) -> dfe_callback_type {
    //            return std::bind(c, std::placeholders::_1);
    //        });

    // py::class_<KTfwd::extensions::gaussian>(m, "GaussianSH")
    //    .def(py::init<double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::gaussian &c) -> dfe_callback_type {
    //            return std::bind(c, std::placeholders::_1);
    //        });

    // py::class_<KTfwd::extensions::gamma>(m, "GammaSH")
    //    .def(py::init<double,double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::gamma &c) -> dfe_callback_type {
    //            return std::bind(c, std::placeholders::_1);
    //        });

    // py::class_<KTfwd::extensions::beta>(m, "BetaSH")
    //    .def(py::init<double,double>())
    //    .def_property_readonly(
    //        "callback",
    //        [](const KTfwd::extensions::beta &c) -> dfe_callback_type {
    //            return std::bind(c, std::placeholders::_1);
    //        });

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
