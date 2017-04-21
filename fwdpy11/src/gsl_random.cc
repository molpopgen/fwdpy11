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
#include <fwdpy11/rng.hpp>
#include <gsl/gsl_randist.h>

namespace py = pybind11;

PYBIND11_PLUGIN(gsl_random)
{
    py::module m("gsl_random",
                 "Wrappers around the GSL random number distributions");

    m.def("gsl_ran_gaussian_ziggurat",
          [](const fwdpy11::GSLrng_t& rng, const double sd) {
              return gsl_ran_gaussian_ziggurat(rng.get(), sd);
          },
          "Gaussian deviate with mean zero");

    m.def("gsl_rng_uniform",
          [](const fwdpy11::GSLrng_t& rng) {
              return gsl_rng_uniform(rng.get());
          },
          "Uniform deviate on interval [0,1)");

    m.def("gsl_ran_flat",
          [](const fwdpy11::GSLrng_t& rng, double a, double b) {
              return gsl_ran_flat(rng.get(), a, b);
          },
          "Unform deviate on interval [a,b)");

    m.def("gsl_ran_poisson",
          [](const fwdpy11::GSLrng_t& rng, double mean) {
              return gsl_ran_poisson(rng.get(), mean);
          },
          "Poisson deviate parameterized by mean.");

    m.def("gsl_ran_geometric",
          [](const fwdpy11::GSLrng_t& rng, const double p) {
              return gsl_ran_geometric(rng.get(), p);
          },
          "Geometric distribution parameterized by success probability.");

    return m.ptr();
}
