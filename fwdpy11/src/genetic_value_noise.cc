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

    virtual pybind11::object
    pickle() const
    {
        return py::make_tuple(mean, sd);
    }

    static inline GaussianNoise
    unpickle(pybind11::object& o)
    {
        py::tuple t(o);
        if (t.size() != 2)
            {
                throw std::runtime_error("invalid object state");
            }
        return GaussianNoise(t[0].cast<double>(), t[1].cast<double>());
    }
};

// Mocking what it takes to get pickling support for a class structure
// and then an object that is composed of unique_ptr to base classes
struct Base
{
    virtual std::string repr() const = 0;
    virtual std::unique_ptr<Base> clone() const = 0;
};

struct Derived : public Base
{
    virtual std::string
    repr() const
    {
        return std::string("fwdpy11.genetic_value_noise.Derived");
    }
    virtual std::unique_ptr<Base>
    clone() const
    {
        return std::unique_ptr<Base>(new Derived(*this));
    }
};

struct Composed
{
    double x;
    std::unique_ptr<Base> b;
    template <typename X> Composed(double a, X&& b_) : x(a), b{ b_.clone() } {}
};

PYBIND11_MODULE(genetic_value_noise, m)
{
    m.doc() = "\"Noise\" added to genetic values.";

    py::class_<fwdpy11::GeneticValueNoise>(
        m, "GeneticValueNoise",
        "ABC for noise classes affecting :class:`fwdpy11.SlocusPop`.");

    py::class_<fwdpy11::NoNoise, fwdpy11::GeneticValueNoise>(
        m, "NoNoise", "Type implying no random effects on genetic values.")
        .def(py::init<>())
        .def(py::pickle(
            [](const fwdpy11::NoNoise& o) -> py::object { return o.pickle(); },
            [](py::object& o) { return fwdpy11::NoNoise::unpickle(o); }));

    py::class_<GaussianNoise, fwdpy11::GeneticValueNoise>(
        m, "GaussianNoise", "Gaussian noise added to genetic values.")
        .def(py::init<double, double>(), py::arg("sd"), py::arg("mean") = 0.0,
             R"delim(
                :param sd: Standard deviation of noise.
                :type sd: float
                :param mean: Mean value of noise.
                :type mean: float
                )delim")
        .def(py::pickle(
            [](const GaussianNoise& o) -> py::object { return o.pickle(); },
            [](py::object& o) { return GaussianNoise::unpickle(o); }));

    // This is an ABC
    py::class_<Base>(m, "Base");

    // This is a concrete class that can be pickled.
    py::class_<Derived, Base>(m, "Derived")
        .def(py::pickle(
            [](const Derived& o) -> py::object { return py::bytes(o.repr()); },
            [](py::object o) -> Derived {
                auto s = o.cast<std::string>();
                if (s.find("Derived") == std::string::npos)
                    {
                        throw std::runtime_error("invalid obect state");
                    }
                return Derived();
            }))
        .def(py::init<>());

    py::class_<Composed>(m, "Composed")
        .def(py::init(
            [](const double d, const Base& b) { return Composed(d, b); }))
        .def_property_readonly("b",
                               [](const Composed& c) { return c.b->repr(); })
        .def(py::pickle(
            // In order to pickle,
            // we return a tuple of c's data plus
            // the result of pickling c.b
            [](const Composed& c) {
                auto x = py::module::import("pickle");
                auto y = x.attr("dumps")(c.b->clone(), -1);
                return py::make_tuple(c.x, y);
            },
            // Unpickling requires some casting magic
            [](py::tuple t) {
                if (t.size() != 2)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto x = py::module::import("pickle");
                auto b = x.attr("loads")(t[1]);
                // We take a const reference to our abstract
                // base in order to compose a new object
                // from the pickled member data.
                // This ONLY works if pickling of classes in
                // our ABC hierarchy is working.
                const Base& br = b.cast<const Base&>();
                return Composed(t[0].cast<double>(), br);
            }));
}
