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
#include <fwdpy11/genetic_values/noise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

#include "wrap_fwdpp_genetic_value.hpp"

namespace py = pybind11;

struct SlocusGBR : public fwdpy11::SlocusPopGeneticValueWithMapping
{
    SlocusGBR(const fwdpy11::GeneticValueIsTrait& g2w)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone() }
    {
    }

    SlocusGBR(const fwdpy11::GeneticValueIsTrait& g2w,
              const fwdpy11::GeneticValueNoise& noise_)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone(),
                                                     noise_.clone() }
    {
    }

    inline double
    sum_haplotype_effect_sizes(const std::size_t gamete_index,
                               const fwdpy11::SlocusPop& pop) const
    {
        double h = 0.0;
        for (auto&& key : pop.gametes[gamete_index].smutations)
            {
                h += pop.mutations[key].s;
            }
        return h;
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        double h1 = sum_haplotype_effect_sizes(
            pop.diploids[diploid_index].first, pop);
        double h2 = sum_haplotype_effect_sizes(
            pop.diploids[diploid_index].second, pop);
        return std::sqrt(h1 * h2);
    }

    inline void
    update(const fwdpy11::SlocusPop& pop)
    {
        gv2w->update(pop);
        noise_fxn->update(pop);
    }
};

// Classes for MlocusPop

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
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, const fwdpy11::GeneticValueIsTrait&>())
        .def(py::init<double, const fwdpy11::GeneticValueIsTrait&,
                      const fwdpy11::GeneticValueNoise&>())
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::SlocusAdditive& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const fwdpy11::SlocusAdditive& wa) {
                return wa.gv.p == fwdpp::additive_diploid::policy::aw;
            });

    py::class_<fwdpy11::SlocusMult, fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusMult")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, const fwdpy11::GeneticValueIsTrait&>())
        .def(py::init<double, const fwdpy11::GeneticValueIsTrait&,
                      const fwdpy11::GeneticValueNoise&>())
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::SlocusMult& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const fwdpy11::SlocusMult& wa) {
                return wa.gv.p == fwdpp::multiplicative_diploid::policy::mw;
            });

    py::class_<SlocusGBR, fwdpy11::SlocusPopGeneticValueWithMapping>(
        m, "SlocusGBR")
        .def(py::init<const fwdpy11::GeneticValueIsTrait&>(), py::arg("gv2w"))
        .def(py::init<const fwdpy11::GeneticValueIsTrait&,
                      const fwdpy11::GeneticValueNoise&>());

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
