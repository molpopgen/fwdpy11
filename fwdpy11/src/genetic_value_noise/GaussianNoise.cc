#include <gsl/gsl_randist.h>
#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

struct GaussianNoise : public fwdpy11::GeneticValueNoise
{
    const double sd, mean;
    GaussianNoise(const double s, const double m) : sd{ s }, mean{ m } {}

    double
    operator()(const fwdpy11::GSLrng_t& rng,
               const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
               const std::size_t /*parent1*/, const std::size_t /*parent2*/,
               const fwdpy11::DiploidPopulation& /*pop*/) const override
    {
        return mean + gsl_ran_gaussian_ziggurat(rng.get(), sd);
    }

    void
    update(const fwdpy11::DiploidPopulation& /*pop*/) override
    {
    }

    std::unique_ptr<fwdpy11::GeneticValueNoise>
    clone() const override
    {
        return std::unique_ptr<GaussianNoise>(new GaussianNoise(*this));
    }
};

void
init_GaussianNoise(py::module& m)
{
    py::class_<GaussianNoise, fwdpy11::GeneticValueNoise>(
        m, "_ll_GaussianNoise")
        .def(py::init<double, double>(), py::arg("sd"), py::arg("mean"));
}

