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
    const double b1, b2, c1, c2;
    // This is our stateful data,
    // which is a record of the
    // additive genetic values of all
    // diploids
    std::vector<double> phenotypes;

    // This constructor is exposed to Python
    snowdrift(double b1_, double b2_, double c1_, double c2_)
        : fwdpy11::DiploidGeneticValue{ 1 }, b1(b1_), b2(b2_), c1(c1_),
          c2(c2_), phenotypes()
    {
    }

    //This constructor makes object unpickling
    //a bit more idiomatic from the C++ point
    //of view.  We implement it as a so-called
    //"perfect-forwarding constructor" to
    //initialize the phenotypes w/o extra copies.
    template <typename T>
    snowdrift(double b1_, double b2_, double c1_, double c2_, T &&p)
        : fwdpy11::DiploidGeneticValue{ 1 }, b1(b1_), b2(b2_), c1(c1_),
          c2(c2_), phenotypes(std::forward<T>(p))
    {
    }

    double
    calculate_gvalue(const std::size_t diploid_index,
                     const fwdpy11::DiploidMetadata& /*metadata*/,
                     const fwdpy11::DiploidPopulation& /*pop*/) const override
    // The call operator must return the genetic value of an individual
    {
        gvalues[0] = phenotypes[diploid_index];
        return gvalues[0];
    }

    double
    genetic_value_to_fitness(
        const fwdpy11::DiploidMetadata &metadata) const override
    // This function converts genetic value to fitness.
    {
        double fitness = 0.0;
        double zself = metadata.g;
        auto N = phenotypes.size();
        for (std::size_t j = 0; j < N; ++j)
            {
                // A record of which diploid we are
                // processesing is the label field of the meta data.
                if (metadata.label != j)
                    {
                        double zpair = zself + phenotypes[j];
                        // Payoff function from Fig 1
                        double a = b1 * zpair + b2 * zpair * zpair - c1 * zself
                                   - c2 * zself * zself;
                        fitness += 1 + std::max(a, 0.0);
                    }
            }
        return fitness / double(N - 1);
    }

    void
    update(const fwdpy11::DiploidPopulation &pop) override
    // A stateful fitness model needs updating.
    {
        phenotypes.resize(pop.N);
        for (std::size_t i = 0; i < pop.N; ++i)
            {
                // A diploid tracks its index via
                // fwdpy11::DiploidMetadata::label
                phenotypes[pop.diploid_metadata[i].label]
                    = fwdpp::additive_diploid(fwdpp::trait(2.0))(
                        pop.diploids[i], pop.haploid_genomes, pop.mutations);
            }
        // This is strictly not necessary in this specific
        // case, but it is required in general, so we
        // do it here by way of example.
        noise_fxn->update(pop);
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
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"))
        .def_readwrite("phenotypes", &snowdrift::phenotypes);
}
