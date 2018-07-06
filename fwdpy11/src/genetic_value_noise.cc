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

#include <gsl/gsl_randist.h>
#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/noise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

struct GaussianNoise : public fwdpy11::GeneticValueNoise
{
    const double sd, mean;
    GaussianNoise(const double s, const double m) : sd{ s }, mean{ m } {}
    virtual double
    operator()(const fwdpy11::GSLrng_t& rng,
               const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
               const std::size_t /*parent1*/, const std::size_t /*parent2*/,
               const fwdpy11::SlocusPop& /*pop*/) const
    {
        return mean + gsl_ran_gaussian_ziggurat(rng.get(), sd);
    }
    virtual double
    operator()(const fwdpy11::GSLrng_t& rng,
               const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
               const std::size_t /*parent1*/, const std::size_t /*parent2*/,
               const fwdpy11::MlocusPop& /*pop*/) const
    {
        return mean + gsl_ran_gaussian_ziggurat(rng.get(), sd);
    }
    DEFAULT_SLOCUSPOP_UPDATE();
    DEFAULT_MLOCUSPOP_UPDATE();
    std::unique_ptr<fwdpy11::GeneticValueNoise>
    clone() const
    {
        return std::unique_ptr<GaussianNoise>(new GaussianNoise(*this));
    }
};

PYBIND11_MODULE(genetic_value_noise, m)
{
    m.doc() = "\"Noise\" added to genetic values.";

    py::class_<fwdpy11::GeneticValueNoise>(
        m, "GeneticValueNoise",
        "ABC for noise classes affecting :class:`fwdpy11.SlocusPop`.");

    py::class_<fwdpy11::NoNoise, fwdpy11::GeneticValueNoise>(
        m, "NoNoise", "Type implying no random effects on genetic values.")
        .def(py::init<>());

    py::class_<GaussianNoise, fwdpy11::GeneticValueNoise>(
        m, "GaussianNoise", "Gaussian noise added to genetic values.")
        .def(py::init<double, double>(), py::arg("sd"), py::arg("mean") = 0.0,
             R"delim(
                :param sd: Standard deviation of noise.
                :type sd: float
                :param mean: Mean value of noise.
                :type mean: float
                )delim");
}
