#include <cmath>
#include <memory>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValueWithMapping.hpp>
#include <fwdpy11/genetic_values/MlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/MlocusPopGeneticValueWithMapping.hpp>
#include <fwdpy11/genetic_values/SlocusAdditive.hpp>
#include <fwdpy11/genetic_values/SlocusMult.hpp>
#include <fwdpy11/genetic_values/SlocusGBR.hpp>
#include <fwdpy11/genetic_values/details/mlocus_aggregators.hpp>
#include <fwdpy11/genetic_values/MlocusAdditive.hpp>
#include <fwdpy11/genetic_values/MlocusMult.hpp>
#include <fwdpy11/genetic_values/MlocusGBR.hpp>
#include <fwdpy11/genetic_values/noise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

PYBIND11_MODULE(genetic_values, m)
{
    auto imported_noise = static_cast<pybind11::object>(
        pybind11::module::import("fwdpy11.genetic_value_noise"));

    py::enum_<fwdpp::additive_diploid::policy>(m, "AdditivePolicy")
        .value("aw", fwdpp::additive_diploid::policy::aw)
        .value("atrait", fwdpp::additive_diploid::policy::atrait)
        .export_values();

    py::enum_<fwdpp::multiplicative_diploid::policy>(m, "MultiplicativePolicy")
        .value("mw", fwdpp::multiplicative_diploid::policy::mw)
        .value("mtrait", fwdpp::multiplicative_diploid::policy::mtrait)
        .export_values();

    py::class_<fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusPopGeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.SlocusPop`");

    py::class_<fwdpy11::SlocusPopGeneticValueWithMapping,
               fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusPopGeneticValueWithMapping")
        .def_property_readonly(
            "gvalue_to_fitness",
            [](const fwdpy11::SlocusPopGeneticValueWithMapping& o) {
                return o.gv2w->clone();
            })
        .def_property_readonly(
            "noise", [](const fwdpy11::SlocusPopGeneticValueWithMapping& o) {
                return o.noise_fxn->clone();
            });

    py::class_<fwdpy11::SlocusAdditive,
               fwdpy11::SlocusPopGeneticValueWithMapping>(m, "SlocusAdditive")
        .def(py::init([](const double scaling) {
                 return fwdpy11::SlocusAdditive(fwdpp::additive_diploid(
                     scaling, fwdpp::additive_diploid::policy::aw));
             }),
             py::arg("scaling"))
        .def(py::init(
            [](const double scaling, const fwdpy11::GeneticValueIsTrait& g) {
                return fwdpy11::SlocusAdditive(
                    fwdpp::additive_diploid(
                        scaling, fwdpp::additive_diploid::policy::atrait),
                    g);
            }))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n) {
            return fwdpy11::SlocusAdditive(
                fwdpp::additive_diploid(
                    scaling, fwdpp::additive_diploid::policy::atrait),
                g, n);
        }))
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::SlocusAdditive& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const fwdpy11::SlocusAdditive& wa) {
                return wa.gv.p == fwdpp::additive_diploid::policy::aw;
            });

    py::class_<fwdpy11::SlocusMult, fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusMult")
        .def(py::init([](const double scaling) {
                 return fwdpy11::SlocusMult(fwdpp::multiplicative_diploid(
                     scaling, fwdpp::multiplicative_diploid::policy::mw));
             }),
             py::arg("scaling"))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g) {
            return fwdpy11::SlocusMult(
                fwdpp::multiplicative_diploid(
                    scaling, fwdpp::multiplicative_diploid::policy::mtrait),
                g);
        }))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n) {
            return fwdpy11::SlocusMult(
                fwdpp::multiplicative_diploid(
                    scaling, fwdpp::multiplicative_diploid::policy::mtrait),
                g, n);
        }))
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::SlocusMult& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const fwdpy11::SlocusMult& wa) {
                return wa.gv.p == fwdpp::multiplicative_diploid::policy::mw;
            });

    py::class_<fwdpy11::SlocusGBR, fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusGBR")
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w) {
            return fwdpy11::SlocusGBR(fwdpy11::GBR(), gv2w);
        }))
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
            return fwdpy11::SlocusGBR(fwdpy11::GBR(), gv2w, noise);
        }));

    // Classes for Mlocus genetic values
    py::class_<fwdpy11::MlocusPopGeneticValue>(
        m, "MlocusPopGeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.MlocusPop`");

    py::class_<fwdpy11::MlocusPopGeneticValueWithMapping,
               fwdpy11::MlocusPopGeneticValue>(
        m, "MlocusPopGeneticValueWithMapping")
        .def_property_readonly(
            "gvalue_to_fitness",
            [](const fwdpy11::MlocusPopGeneticValueWithMapping& o) {
                return o.gv2w->clone();
            })
        .def_property_readonly(
            "noise", [](const fwdpy11::MlocusPopGeneticValueWithMapping& o) {
                return o.noise_fxn->clone();
            });

    py::class_<fwdpy11::MlocusAdditive>(m, "MlocusAdditive")
        .def(py::init([](const double scaling) {
                 return fwdpy11::MlocusAdditive(
                     fwdpp::additive_diploid(
                         scaling, fwdpp::additive_diploid::policy::aw),
                     fwdpy11::aggregate_additive_fitness());
             }),
             py::arg("scaling"))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& gv2w) {
            return fwdpy11::MlocusAdditive(
                fwdpp::additive_diploid(
                    scaling, fwdpp::additive_diploid::policy::atrait),
                fwdpy11::aggregate_additive_fitness(), gv2w);
        }))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
            return fwdpy11::MlocusAdditive(
                fwdpp::additive_diploid(
                    scaling, fwdpp::additive_diploid::policy::atrait),
                fwdpy11::aggregate_additive_fitness(), gv2w, noise);
        }));

    py::class_<fwdpy11::MlocusMult>(m, "MlocusMult")
        .def(py::init([](const double scaling) {
                 return fwdpy11::MlocusMult(
                     fwdpp::multiplicative_diploid(
                         scaling, fwdpp::multiplicative_diploid::policy::mw),
                     fwdpy11::aggregate_mult_fitness());
             }),
             py::arg("scaling"))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& gv2w) {
            return fwdpy11::MlocusMult(
                fwdpp::multiplicative_diploid(
                    scaling, fwdpp::multiplicative_diploid::policy::mtrait),
                fwdpy11::aggregate_mult_fitness(), gv2w);
        }))
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
            return fwdpy11::MlocusMult(
                fwdpp::multiplicative_diploid(
                    scaling, fwdpp::multiplicative_diploid::policy::mtrait),
                fwdpy11::aggregate_mult_fitness(), gv2w, noise);
        }));

    py::class_<fwdpy11::MlocusGBR>(m, "MlocusGBR")
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w) {
            return fwdpy11::MlocusGBR(
                fwdpy11::GBR(), fwdpy11::aggregate_additive_trait(), gv2w);
        }))
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
            return fwdpy11::MlocusGBR(fwdpy11::GBR(),
                                      fwdpy11::aggregate_additive_trait(),
                                      gv2w, noise);
        }));

    py::class_<fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueToFitnessMap",
        "ABC for functions translating genetic values into fitness.");

    py::class_<fwdpy11::GeneticValueIsTrait,
               fwdpy11::GeneticValueToFitnessMap>(m, "GeneticValueIsTrait",
                                                  "ABC");

    py::class_<fwdpy11::GeneticValueIsFitness,
               fwdpy11::GeneticValueToFitnessMap>(m, "GeneticValueIsFitness")
        .def(py::init<>());

    // TODO: need to decide on (VS, opt) vs (opt, VS) and have
    // GSS and GSSmo do the same thing
    py::class_<fwdpy11::GSS, fwdpy11::GeneticValueIsTrait>(
        m, "GSS", "Gaussian stabilizing selection.")
        .def(py::init<double, double>(), py::arg("VS"), py::arg("opt"))
        .def_readonly("VS", &fwdpy11::GSS::VS)
        .def_readonly("opt", &fwdpy11::GSS::opt);

    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "GSSmo", "Gaussian stabilizing selection with a moving optimum.")
        .def(
            py::init<std::vector<std::tuple<std::uint32_t, double, double>>>(),
            py::arg("optima"));
}
