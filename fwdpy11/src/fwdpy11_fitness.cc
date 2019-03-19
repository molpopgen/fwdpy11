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
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/fitness/single_locus_stateless_fitness.hpp>

namespace py = pybind11;

PYBIND11_MODULE(fitness, m)
{
    m.doc() = "Fitness models.";

    py::class_<fwdpy11::single_locus_fitness,
               std::shared_ptr<fwdpy11::single_locus_fitness>>(m,
                                                               "DiploidFitness",
                                                               R"delim(
            A fitness function or trait value function
            for a single-deme, single-region simulation (
            :class:`fwdpy11.Diploid`).
            
            Current fitness ovjects derived from this type are:
            
            * :class:`fwdpy11.fitness.DiploidAdditive`
            * :class:`fwdpy11.fitness.DiploidMult`
            
            Current trait value objects derived from this type are:

            * :class:`fwdpy11.trait_values.DiploidAdditiveTrait`
            * :class:`fwdpy11.trait_values.DiploidMultTrait`
            * :class:`fwdpy11.trait_values.DiploidGBRTrait`
            )delim")
        .def("__call__",
             [](const std::shared_ptr<fwdpy11::single_locus_fitness>& aw,
                const fwdpy11::Diploid& dip,
                const fwdpy11::DiploidPopulation& pop) {
                 return aw->callback()(dip, pop.gametes, pop.mutations);
             });

    // pybind11::class_<fwdpy11::single_locus_stateless_fitness,
    //                 std::shared_ptr<fwdpy11::single_locus_stateless_fitness>,
    //                 fwdpy11::single_locus_fitness>(
    //    m, "DiploidCustomStatelessGeneticValue", "Custom stateless genetic "
    //                                            "value/fitness function. See
    //                                            "
    //                                            ":ref:`customgvaluecpp`.")
    //    .def(py::init<fwdpy11::single_locus_fitness_fxn>());

    py::class_<fwdpy11::single_locus_mult_wrapper,
               std::shared_ptr<fwdpy11::single_locus_mult_wrapper>,
               fwdpy11::single_locus_fitness>(m, "DiploidMult", R"delim(
        Multiplicative fitness for single-deme simulations.
        At a single mutation, fitness is 0, 1+sh, 1+scaling*s
        for genotypes AA, Aa, and aa, respectively. The scaling
        parameter is a constructor argument.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.DiploidMult(1.0)
                       )delim")
        .def(py::init<double>(), py::arg("scaling") = 2.0)
        .def_readonly("scaling", &fwdpy11::single_locus_mult_wrapper::scaling,
                      "Get the scaling parameter.")
        .def("__repr__",
             [](const fwdpy11::single_locus_mult_wrapper& m) {
                 std::string rv = "fitness.DiploidMult(";
                 rv += std::to_string(m.scaling);
                 rv += ')';
                 return rv;
             })
        .def(py::pickle(
            [](const fwdpy11::single_locus_mult_wrapper& w) {
                return py::make_tuple(w.scaling);
            },
            [](py::tuple t) {
                return std::make_shared<fwdpy11::single_locus_mult_wrapper>(
                    t[0].cast<double>());
            }));

    py::class_<fwdpy11::single_locus_additive_wrapper,
               std::shared_ptr<fwdpy11::single_locus_additive_wrapper>,
               fwdpy11::single_locus_fitness>(m, "DiploidAdditive", R"delim(
        Additive fitness for single-deme simulations.
        Fitness is max(0,1 + :math:`\sum_{i} x_i`), 
        where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
        for AA, Aa, and aa, respectively.

        .. testcode::

            import fwdpy11.fitness as fp11w
            w = fp11w.DiploidAdditive(2.0)
        )delim")
        .def_readonly("scaling",
                      &fwdpy11::single_locus_additive_wrapper::scaling,
                      "Get the scaling parameter.")
        .def(py::init<double>(), py::arg("scaling") = 2.0)
        .def("__repr__",
             [](const fwdpy11::single_locus_additive_wrapper& m) {
                 std::string rv = "fitness.DiploidAdditive(";
                 rv += std::to_string(m.scaling);
                 rv += ')';
                 return rv;
             })
        .def(py::pickle(
            [](const fwdpy11::single_locus_additive_wrapper& w) {
                return py::make_tuple(w.scaling);
            },
            [](py::tuple t) {
                return std::
                    make_shared<fwdpy11::single_locus_additive_wrapper>(
                        t[0].cast<double>());
            }));
}
