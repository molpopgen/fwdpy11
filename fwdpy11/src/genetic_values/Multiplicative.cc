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

    using DiploidMultiplicative
        = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
            single_deme_multiplicative_het, single_deme_multiplicative_hom,
            multi_deme_multiplicative_het, multi_deme_multiplicative_hom, 1>;
}

void
init_Multiplicative(py::module& m)
{
    py::class_<DiploidMultiplicative, fwdpy11::DiploidGeneticValue>(m,
                                                                    "_ll_Multiplicative")
        .def(py::init([](double scaling,
                         const fwdpy11::GeneticValueIsTrait* gvalue_to_fitness,
                         const fwdpy11::GeneticValueNoise* noise, std::size_t ndemes) {
                 if (gvalue_to_fitness != nullptr)
                     {
                         return DiploidMultiplicative(ndemes, scaling,
                                                      final_multiplicative_trait(),
                                                      gvalue_to_fitness, noise);
                     }
                 return DiploidMultiplicative(ndemes, scaling,
                                              final_multiplicative_fitness(),
                                              gvalue_to_fitness, noise);
             }),
             py::arg("scaling"), py::arg("gvalue_to_fitness"), py::arg("noise"),
             py::arg("ndemes"));
}
