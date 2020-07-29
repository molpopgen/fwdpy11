/* Implement a stateful fitness model.
 * We define a new C++ type that will be
 * wrapped as a fwdpy11.DiploidGeneticValue
 * object.
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>

namespace py = pybind11;

struct snowdrift : public fwdpy11::DiploidGeneticValue
/* This is the low-level implementation of our
 * stateful fitness object.  It records the
 * model parameters and holds a
 * vector to track individual phenotypes.
 *
 * Here, we publicly inherit from fwdpy11::DiploidGeneticValue,
 * which is defined in the header included above.  It is
 * an abstract class in C++ terms, and is reflected
 * as a Python Abstract Base Class (ABC) called
 * fwdpy11.DiploidGeneticValue.
 *
 * The phenotypes get updated each generation during
 * the simulation.
 *
 * The phenotypes will be the simple additive model,
 * calculated using fwdpp's machinery.
 */
{
    const double b1, b2, c1, c2, slope, sig0;
    // This is our stateful data,
    // which is a record of the
    // additive genetic values of all
    // diploids
    std::vector<double> phenotypes;
    const fwdpp::additive_diploid additive;

    // This constructor is exposed to Python
    snowdrift(double b1_, double b2_, double c1_, double c2_, double slope, double p0)
        : fwdpy11::DiploidGeneticValue{1, nullptr, nullptr}, b1(b1_), b2(b2_), c1(c1_),
          c2(c2_), slope(slope), sig0(1. / slope * std::log(p0 / (1. - p0))),
          phenotypes(), additive{fwdpp::trait{2.0}}
    {
        if (p0 == 0 || p0 == 1.)
            {
                throw std::invalid_argument("p0 must be 0 < p0 < 1");
            }
    }

    double
    calculate_gvalue(const fwdpy11::DiploidGeneticValueData data) override
    // The call operator must return the genetic value of an individual
    {
        gvalues[0] = phenotypes[data.metadata_index];
        return gvalues[0];
    }

    double
    genetic_value_to_fitness(
        const fwdpy11::DiploidGeneticValueToFitnessData data) override
    // This function converts genetic value to fitness.
    {
        double zself = data.offspring_metadata.get().g;
        auto N = phenotypes.size();
        auto other = gsl_rng_uniform_int(data.rng.get().get(), N);
        while (other == data.offspring_metadata.get().label)
            {
                other = gsl_rng_uniform_int(data.rng.get().get(), N);
            }
        double zpair = zself + phenotypes[other];
        double a
            = 1. + b1 * zpair + b2 * zpair * zpair - c1 * zself - c2 * zself * zself;
        return std::max(a, 0.0);
    }

    void
    update(const fwdpy11::DiploidPopulation &pop) override
    // A stateful fitness model needs updating.
    {
        phenotypes.clear();
        for (auto &md : pop.diploid_metadata)
            {
                // A diploid tracks its index via
                // fwdpy11::DiploidMetadata::label
                auto g = additive(pop.diploids[md.label], pop.haploid_genomes,
                                  pop.mutations);
                phenotypes.push_back(1. / (1. + std::exp(-slope * (g + sig0))));
            }
    }
};

PYBIND11_MODULE(ll_snowdrift, m)
{
    m.doc() = "Example of custom stateful fitness model.";

    // We need to import the Python version of our base class:
    pybind11::object imported_snowdrift_base_class_type
        = pybind11::module::import("fwdpy11").attr("DiploidGeneticValue");

    // Create a Python class based on our new type
    py::class_<snowdrift, fwdpy11::DiploidGeneticValue>(m, "_ll_DiploidSnowdrift")
        .def(py::init<double, double, double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"), py::arg("slope"),
             py::arg("p0"))
        .def_readwrite("phenotypes", &snowdrift::phenotypes);
}
