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
#include <functional>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/genetic_values/fwdpp_wrappers/fwdpp_genetic_value.hpp>
#include <fwdpy11/types/Mutation.hpp>
#include <pybind11/pybind11.h>
#include "gvalue_pickle_helpers.hpp"

namespace py = pybind11;

namespace
{
    struct single_deme_multiplicative_het
    {
        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + m.s * m.h);
        }
    };

    struct multi_deme_multiplicative_het
    {
        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + m.esizes[deme] * m.heffects[deme]);
        }
    };

    struct single_deme_multiplicative_hom
    {
        double scaling;
        single_deme_multiplicative_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + scaling * m.s);
        }
    };

    struct multi_deme_multiplicative_hom
    {
        double scaling;
        multi_deme_multiplicative_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d *= (1. + scaling * m.esizes[deme]);
        }
    };

    struct final_multiplicative_trait
    {
        inline double
        operator()(double d) const
        {
            return d - 1.0;
        }
    };

    struct final_multiplicative_fitness
    {
        inline double
        operator()(double d) const
        {
            return std::max(0.0, d);
        }
    };

    using DiploidMult = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
        single_deme_multiplicative_het, single_deme_multiplicative_hom,
        multi_deme_multiplicative_het, multi_deme_multiplicative_hom, pickle_gvalue, 1>;
} // namespace fwdpy11

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
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
)delim";

static const auto MULT_CONSTRUCTOR_3 =
    R"delim(
Multiplicative effects on a trait with a specific mapping from 
genetic value to fitness and random effects ("noise").

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
:param noise: Function to generate random effects on trait value.
:type noise: :class:`fwdpy11.GeneticValueNoise`
)delim";

void
init_Multiplicative(py::module& m)
{
    py::class_<DiploidMult, fwdpy11::DiploidGeneticValue>(
        m, "Multiplicative", "Multiplicative genetic values.")
        .def(py::init([](const double scaling, std::size_t ndemes) {
                 return DiploidMult(ndemes, scaling, final_multiplicative_fitness());
             }),
             py::arg("scaling"), py::arg("ndemes") = 1, MULT_CONSTRUCTOR_1)
        .def(py::init([](const double scaling, const fwdpy11::GeneticValueIsTrait& g,
                         std::size_t ndemes) {
                 return DiploidMult(ndemes, scaling, final_multiplicative_trait(), g);
             }),
             py::arg("scaling"), py::arg("gv2w"), py::arg("ndemes") = 1,
             MULT_CONSTRUCTOR_2)
        .def(py::init([](const double scaling, const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n, std::size_t ndemes) {
                 return DiploidMult(ndemes, scaling, final_multiplicative_trait(), g, n);
             }),
             py::arg("scaling"), py::arg("gv2w"), py::arg("noise"),
             py::arg("ndemes") = 1, MULT_CONSTRUCTOR_3)
        .def_property_readonly("scaling", &DiploidMult::scaling,
                               "Access to the scaling parameter.")
        .def_property_readonly(
            "is_fitness",
            [](const DiploidMult& self) {
                PyErr_WarnEx(
                    PyExc_DeprecationWarning,
                    "Additive.is_fitness is deprecated.  Use maps_to_fitness instead",
                    0);
                return self.is_fitness();
            },
            R"delim(
            Returns True if instance calculates fitness as the genetic value
            and False if the genetic value is a trait value.
            
            .. deprecated:: 0.7.0

                Use :attr:`fwdpy11.DiploidGeneticValue.maps_to_fitness` instead.
            )delim")

        .def(py::pickle(
            [](const DiploidMult& a) {
                auto p = py::module::import("pickle");
                return py::make_tuple(a.pickle(), p.attr("dumps")(a.gv2w->clone(), -1),
                                      p.attr("dumps")(a.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                return unpickle_gvalue<DiploidMult, final_multiplicative_fitness,
                                       final_multiplicative_trait, std::true_type>(t);
            }));
}
