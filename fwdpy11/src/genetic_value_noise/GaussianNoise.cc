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
               const fwdpy11::DiploidPopulation& /*pop*/) const
    {
        return mean + gsl_ran_gaussian_ziggurat(rng.get(), sd);
    }

    DEFAULT_DIPLOID_POP_UPDATE();
    
    std::unique_ptr<fwdpy11::GeneticValueNoise>
    clone() const
    {
        return std::unique_ptr<GaussianNoise>(new GaussianNoise(*this));
    }

    virtual pybind11::object
    pickle() const
    {
        return py::make_tuple(sd, mean);
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

void
init_GaussianNoise(py::module& m)
{
    py::class_<GaussianNoise, fwdpy11::GeneticValueNoise>(
        m, "GaussianNoise", "Gaussian noise added to genetic values.")
        .def(py::init<double, double>(), py::arg("sd"), py::arg("mean") = 0.0,
             R"delim(
                :param sd: Standard deviation of noise.
                :type sd: float
                :param mean: Mean value of noise.
                :type mean: float
                )delim")
        .def_readonly("mean", &GaussianNoise::mean)
        .def_readonly("sd", &GaussianNoise::sd)
        .def(py::pickle(
            [](const GaussianNoise& o) -> py::object { return o.pickle(); },
            [](py::object& o) { return GaussianNoise::unpickle(o); }));
}

