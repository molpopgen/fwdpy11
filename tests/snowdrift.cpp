/* Implement a stateful fitness model.
 * We define a new C++ type that will be
 * wrapped as a fwdpy11.fitness.DiploidFitness
 * object.
 *
 * Such a fitness model is ultimately responsible
 * for generating a bound C++ callback whose signature
 * is fwdpy11::single_locus_fitness_fxn.
 *
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>

namespace py = pybind11;

struct snowdrift : public fwdpy11::DiploidPopulationGeneticValue
/* This is our stateful fitness object.
 * It records the model parameters and holds a
 * vector to track individual phenotypes.
 *
 * Here, we publicly inherit from fwdpy11::DiploidPopulationGeneticValue,
 * which is defined in the header included above.  It is
 * an abstract class in C++ terms, and is reflected
 * as a Python Abstract Base Class (ABC) called
 * fwdpy11.genetic_values.DiploidPopulationGeneticValue.
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
        : fwdpy11::DiploidPopulationGeneticValue{ 1 }, b1(b1_), b2(b2_),
          c1(c1_), c2(c2_), phenotypes()
    {
    }

    //This constructor makes object unpickling
    //a bit more idiomatic from the C++ point
    //of view.  We implement it as a so-called
    //"perfect-forwarding constructor" to
    //initialize the phenotypes w/o extra copies.
    template <typename T>
    snowdrift(double b1_, double b2_, double c1_, double c2_, T &&p)
        : fwdpy11::DiploidPopulationGeneticValue{ 1 }, b1(b1_), b2(b2_),
          c1(c1_), c2(c2_), phenotypes(std::forward<T>(p))
    {
    }

    inline double
    calculate_gvalue(const std::size_t diploid_index,
                     const fwdpy11::DiploidPopulation & /*pop*/) const
    // The call operator must return the genetic value of an individual
    {
        gvalues[0] = phenotypes[diploid_index];
        return gvalues[0];
    }

    inline double
    genetic_value_to_fitness(const fwdpy11::DiploidMetadata &metadata) const
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

    inline double
    noise(const fwdpy11::GSLrng_t & /*rng*/,
          const fwdpy11::DiploidMetadata & /*offspring_metadata*/,
          const std::size_t /*parent1*/, const std::size_t /*parent2*/,
          const fwdpy11::DiploidPopulation & /*pop*/) const
    // This function may be used to model random effects...
    {
        //...but there are no random effects here.
        return 0.0;
    }

    inline void
    update(const fwdpy11::DiploidPopulation &pop)
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
    }

    // In order to support pickling, the ABC requries
    // that the following function be defined for a subclass.
    // It must return the relevant data needed for serialization
    // as a Python object.  pybind11 makes this easy (for most
    // cases, most of the time).
    py::object
    pickle() const
    {
        return py::make_tuple(b1, b2, c1, c2, phenotypes);
    }

    py::tuple
    shape() const
    {
        return py::make_tuple(1);
    }
};

PYBIND11_MODULE(snowdrift, m)
{
    m.doc() = "Example of custom stateful fitness model.";

    // We need to import the Python version of our base class:
    pybind11::object imported_snowdrift_base_class_type
        = pybind11::module::import("fwdpy11").attr("GeneticValue");

    // Create a Python class based on our new type
    py::class_<snowdrift, fwdpy11::DiploidPopulationGeneticValue>(
        m, "DiploidSnowdrift")
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"))
        .def_readonly("b1", &snowdrift::b1)
        .def_readonly("b2", &snowdrift::b2)
        .def_readonly("c1", &snowdrift::c1)
        .def_readonly("c2", &snowdrift::c2)
        .def_readwrite("phenotypes", &snowdrift::phenotypes)
        // Implement pickling support
        .def(py::pickle(
            //Pickling here is quite simple.  We simply
            //return the results of the member fxn
            [](const snowdrift &s) { return s.pickle(); },
            //Unpickling is almost always harder:
            //Note that any of the Python steps below
            //will raise exceptions if they fail,
            //making this code safe at run time.
            [](py::object o) {
                //Convert object to tuple.
                py::tuple t(o);
                // Check tuple size and
                // throw error if it is not right
                if (t.size() != 5)
                    {
                        throw std::runtime_error("invalid object state");
                    }

                //Get our data out via pybind11's casting
                //mechanisms.  If things are not
                //cast-able, exceptions are thrown
                double b1 = t[0].cast<double>();
                double b2 = t[1].cast<double>();
                double c1 = t[2].cast<double>();
                double c2 = t[3].cast<double>();
                auto p = t[4].cast<std::vector<double>>();
                //Create our object, using move semantics
                //for efficiency
                snowdrift rv(b1, b2, c1, c2, std::move(p));
                return rv;
            }));
}
