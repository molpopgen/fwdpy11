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
#include <fwdpy11/genetic_values/noise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

// Define some docstrings as static variable to avoid
// repetition below

void init_genetic_values(py::module&);

PYBIND11_MODULE(genetic_values, m)
{
    auto imported_noise = static_cast<pybind11::object>(
        pybind11::module::import("fwdpy11.genetic_value_noise"));
    init_genetic_values(m);

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

}
