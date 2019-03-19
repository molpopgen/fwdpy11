//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <cmath>
#include <memory>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValueWithMapping.hpp>
#include <fwdpy11/genetic_values/DiploidAdditive.hpp>
#include <fwdpy11/genetic_values/DiploidMult.hpp>
#include <fwdpy11/genetic_values/DiploidGBR.hpp>
#include <fwdpy11/genetic_values/details/mlocus_aggregators.hpp>
#include <fwdpy11/genetic_values/noise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

// Define some docstrings as static variable to avoid
// repetition below

static const auto GBR_CONSTRUCTOR1 =
    R"delim(
 Construct object with specific genetic value to fitness map.
 
 :param gv2w: Genetic value to fitness map
 :type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
 )delim";

static const auto GBR_CONSTRUCTOR2 =
    R"delim(
Construct object with specific genetic value to fitness map 
and random effects on trait value.

:param gv2w: Genetic value to fitness map
:type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
:param noise: Model of random effects on trait value.
:type noise: :class:`fwdpy11.genetic_value_noise.GeneticValueNoise`
)delim";

static const auto ADDITIVE_CONSTRUCTOR_1 =
    R"delim(
Additive effects on fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float

For a model of fitness, the genetic value is 1, 1+e*h,
1+scaling*e for genotypes AA, Aa, and aa, respectively.
)delim";

static const auto ADDITIVE_CONSTRUCTOR_2 =
    R"delim(
Construct an object of additive effects on a trait with a specific
functional mapping from genetic value to fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
)delim";

static const auto ADDITIVE_CONSTRUCTOR_3 =
    R"delim(
Additive effects on a trait with a specific mapping from 
genetic value to fitness and random effects ("noise").

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
:param noise: Function to generate random effects on trait value.
:type noise: :class:`fwdpy11.genetic_value_noise.GeneticValueNoise`
)delim";

static const auto MULT_CONSTRUCTOR_1 =
    R"delim(
Multiplicative effects on fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float

For a model of fitness, the genetic value is 1, 1+e*h,
1+scaling*e for genotypes AA, Aa, and aa, respectively.
)delim";

static const auto MULT_CONSTRUCTOR_2 =
    R"delim(
Construct an object of multiplicative effects on a trait with a specific
functional mapping from genetic value to fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
)delim";

static const auto MULT_CONSTRUCTOR_3 =
    R"delim(
Multiplicative effects on a trait with a specific mapping from 
genetic value to fitness and random effects ("noise").

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.genetic_values.GeneticValueIsTrait`
:param noise: Function to generate random effects on trait value.
:type noise: :class:`fwdpy11.genetic_value_noise.GeneticValueNoise`
)delim";

void init_genetic_values(py::module&);

PYBIND11_MODULE(genetic_values, m)
{
    auto imported_noise = static_cast<pybind11::object>(
        pybind11::module::import("fwdpy11.genetic_value_noise"));

    py::class_<fwdpy11::DiploidPopulationGeneticValue>(
        m, "DiploidPopulationGeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.DiploidPopulation`")
        .def("__call__",
             [](const fwdpy11::DiploidPopulationGeneticValue& gv,
                const std::size_t diploid_index,
                const fwdpy11::DiploidPopulation& pop) {
                 return gv.calculate_gvalue(diploid_index, pop);
             },
             R"delim(
             :param diploid_index: The index of the individual to calculate.
             :type diploid_index: int >= 0
             :param pop: The population object containing the individual.
             :type pop: :class:`fwdpy11.DiploidPopulation`
             :return: The genetic value of an individual.
             :rtype: float
             )delim",
             py::arg("diploid_index"), py::arg("pop"))
        .def("fitness",
             [](const fwdpy11::DiploidPopulationGeneticValue& gv,
                const std::size_t diploid_index,
                const fwdpy11::DiploidPopulation& pop) {
                 return gv.genetic_value_to_fitness(
                     pop.diploid_metadata[diploid_index]);
             },
             R"delim(
        :param diploid_index: The index of the individual
        :type diploid_index: int >= 0
        :param pop: The population containing the individual
        :type pop: :class:`fwdpy11.DiploidPopulation`
        :return: The fitness of an individual.
        :rtype: float
        )delim",
             py::arg("diploid_index"), py::arg("pop"))
        .def_property_readonly("shape",
                               [](const fwdpy11::DiploidPopulationGeneticValue& self) {
                                   return self.shape();
                               },
                               R"delim(
                               Return the dimensions of the genetic values.

                               .. versionadded:: 0.3.0
                               )delim")
        .def_readonly("genetic_values",
                      &fwdpy11::DiploidPopulationGeneticValue::gvalues,
                      R"delim(
                      Return the list of genetic values.

                      .. versionadded:: 0.3.0
                      )delim");

    py::class_<fwdpy11::DiploidPopulationGeneticValueWithMapping,
               fwdpy11::DiploidPopulationGeneticValue>(
        m, "DiploidPopulationGeneticValueWithMapping",
        "Genetic value calculations with flexible mapping of genetic value to "
        "fitness.")
        .def_property_readonly(
            "gvalue_to_fitness",
            [](const fwdpy11::DiploidPopulationGeneticValueWithMapping& o) {
                return o.gv2w->clone();
            },
            "Access the genetic value to fitness map.")
        .def_property_readonly(
            "noise",
            [](const fwdpy11::DiploidPopulationGeneticValueWithMapping& o) {
                return o.noise_fxn->clone();
            },
            "Access the random noise funcion");

    py::class_<fwdpy11::DiploidAdditive,
               fwdpy11::DiploidPopulationGeneticValueWithMapping>(
        m, "DiploidAdditive", "Additive genetic values.")
        .def(py::init([](const double scaling) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::fitness(scaling)));
             }),
             py::arg("scaling"), ADDITIVE_CONSTRUCTOR_1)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::trait(scaling)), g);
             }),
             py::arg("scaling"), py::arg("gv2w"), ADDITIVE_CONSTRUCTOR_2)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::trait(scaling)), g, n);
             }),
             py::arg("scaling"), py::arg("gv2w"), py::arg("noise"),
             ADDITIVE_CONSTRUCTOR_3)
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::DiploidAdditive& wa) { return wa.gv.scaling; },
            "Access to the scaling parameter.")
        .def_property_readonly(
            "is_fitness",
            [](const fwdpy11::DiploidAdditive& wa) {
                return wa.gv.gvalue_is_fitness;
            },
            "Returns True if instance calculates fitness as the genetic value "
            "and False if the genetic value is a trait value.")
        .def(py::pickle(
            [](const fwdpy11::DiploidAdditive& a) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    a.pickle(), p.attr("dumps")(a.gv2w->clone(), -1),
                    p.attr("dumps")(a.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto t0 = t[0].cast<py::tuple>();
                int pol = t0[0].cast<int>();

                double scaling = t0[1].cast<double>();
                auto p = py::module::import("pickle");
                auto a
                    = (pol == 1)
                          ? fwdpp::additive_diploid(fwdpp::trait(scaling))
                          : fwdpp::additive_diploid(fwdpp::fitness(scaling));
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::DiploidAdditive(
                    std::move(a),
                    t1.cast<const fwdpy11::GeneticValueToFitnessMap&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));

    py::class_<fwdpy11::DiploidMult, fwdpy11::DiploidPopulationGeneticValueWithMapping>(
        m, "DiploidMult", "Multiplicative genetic values.")
        .def(py::init([](const double scaling) {
                 return fwdpy11::DiploidMult(
                     fwdpp::multiplicative_diploid(fwdpp::fitness(scaling)));
             }),
             py::arg("scaling"), MULT_CONSTRUCTOR_1)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g) {
                 return fwdpy11::DiploidMult(
                     fwdpp::multiplicative_diploid(fwdpp::trait(scaling)), g);
             }),
             py::arg("scaling"), py::arg("gv2w"), MULT_CONSTRUCTOR_2)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n) {
                 return fwdpy11::DiploidMult(
                     fwdpp::multiplicative_diploid(fwdpp::trait(scaling)), g,
                     n);
             }),
             py::arg("scaling"), py::arg("gv2w"), py::arg("noise"),
             MULT_CONSTRUCTOR_3)
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::DiploidMult& wa) { return wa.gv.scaling; },
            "Access to the scaling parameter.")
        .def_property_readonly(
            "is_fitness",
            [](const fwdpy11::DiploidMult& wa) {
                return wa.gv.gvalue_is_fitness;
            },
            "Returns True if instance calculates fitness as the genetic "
            "value "
            "and False if the genetic value is a trait value.")
        .def(py::pickle(
            [](const fwdpy11::DiploidMult& a) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    a.pickle(), p.attr("dumps")(a.gv2w->clone(), -1),
                    p.attr("dumps")(a.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto t0 = t[0].cast<py::tuple>();
                int pol = t0[0].cast<int>();

                double scaling = t0[1].cast<double>();
                auto p = py::module::import("pickle");
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                auto a = (pol == 1) ? fwdpp::multiplicative_diploid(
                                          fwdpp::trait(scaling))
                                    : fwdpp::multiplicative_diploid(
                                          fwdpp::fitness(scaling));
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::DiploidMult(
                    std::move(a),
                    t1.cast<const fwdpy11::GeneticValueToFitnessMap&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));

    py::class_<fwdpy11::DiploidGBR, fwdpy11::DiploidPopulationGeneticValueWithMapping>(
        m, "DiploidGBR",
        R"delim(
        The "gene-based recessive" trait model described in Thornton et al.
        2013 http://dx.doi.org/10.1371/journal.pgen.1003258 and Sanjak et al. 2017
        http://dx.doi.org/10.1371/journal.pgen.1006573.

        The trait value is the geometric mean of the sum of effect sizes on each haplotype.
        It is undefined for the case where these sums are negative.
        )delim")
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w) {
                 return fwdpy11::DiploidGBR(fwdpy11::GBR{}, gv2w);
             }),
             py::arg("gv2w"), GBR_CONSTRUCTOR1)
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
                 return fwdpy11::DiploidGBR(fwdpy11::GBR{}, gv2w, noise);
             }),
             py::arg("gv2w"), py::arg("noise"), GBR_CONSTRUCTOR2)
        .def(py::pickle(
            [](const fwdpy11::DiploidGBR& g) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    g.pickle(), p.attr("dumps")(g.gv2w->clone(), -1),
                    p.attr("dumps")(g.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                std::string s = t[0].cast<std::string>();
                if (s != "GBR")
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto p = py::module::import("pickle");
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::DiploidGBR(
                    fwdpy11::GBR{},
                    t1.cast<const fwdpy11::GeneticValueIsTrait&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));

    py::class_<fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueToFitnessMap",
        "ABC for functions translating genetic values into fitness.");

    py::class_<fwdpy11::GeneticValueIsTrait,
               fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueIsTrait",
        "ABC for functions mapping genetic values representing traits to "
        "fitness.");

    py::class_<fwdpy11::GeneticValueIsFitness,
               fwdpy11::GeneticValueToFitnessMap>(
        m, "GeneticValueIsFitness",
        "Type implying the the genetic value is fitness.")
        .def(py::init<>())
        .def(py::pickle(
            [](const fwdpy11::GeneticValueIsFitness& g) { return g.pickle(); },
            [](py::object o) {
                std::string s = o.cast<std::string>();
                if (s.find("GeneticValueIsFitness") == std::string::npos)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                return fwdpy11::GeneticValueIsFitness();
            }));

    py::class_<fwdpy11::GSS, fwdpy11::GeneticValueIsTrait>(
        m, "GSS", "Gaussian stabilizing selection.")
        .def(py::init<double, double>(), py::arg("opt"), py::arg("VS"),
             R"delim(
                :param opt: Optimal trait value.
                :type opt: float
                :param VS: Strength of stabilizing selection
                :type VS: float
                )delim")
        .def_readonly("VS", &fwdpy11::GSS::VS, "Read-only access to VS")
        .def_readonly("opt", &fwdpy11::GSS::opt,
                      "Read-only access to optimal trait value.")
        .def(py::pickle(
            [](const fwdpy11::GSS& g) { return g.pickle(); },
            [](py::object o) {
                py::tuple t(o);
                if (t.size() != 2)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                return fwdpy11::GSS(t[0].cast<double>(), t[1].cast<double>());
            }));

    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "GSSmo", "Gaussian stabilizing selection with a moving optimum.")
        .def(
            py::init<std::vector<std::tuple<std::uint32_t, double, double>>>(),
            py::arg("optima"),
            R"delim(
            :param optima: Model parameters over time
            :type optima: list
            
            Each element of optima must be a tuple of 
            (generation, optimal trait value, VS)
            )delim")
        .def_readonly("VS", &fwdpy11::GSSmo::VS)
        .def_readonly("opt", &fwdpy11::GSSmo::opt)
        .def_readonly("optima", &fwdpy11::GSSmo::optima)
        .def(py::pickle(
            [](const fwdpy11::GSSmo& g) { return g.pickle(); },
            [](py::object o) {
                py::tuple t(o);
                if (t.size() != 4)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto opt = t[0].cast<double>();
                auto VS = t[1].cast<double>();
                auto co
                    = t[2].cast<decltype(fwdpy11::GSSmo::current_optimum)>();
                auto optima = t[3].cast<decltype(fwdpy11::GSSmo::optima)>();
                auto rv = fwdpy11::GSSmo(std::move(optima));
                rv.opt = opt;
                rv.VS = VS;
                rv.current_optimum = co;
                return rv;
            }));

    init_genetic_values(m);
}
