/* Implement a stateful fitness model.
 * We define a new C++ type that will be
 * wrapped as a fwdpy11.fitness.SlocusFitness
 * object.
 *
 * Such a fitness model is ultimately responsible
 * for generating a bound C++ callback whose signature
 * is fwdpy11::single_locus_fitness_fxn.
 *
 * The module is built using cppimport:
 * https://github.com/tbenthompson/cppimport
 */

/* The next block of code is used by cppimport
 * The formatting is important, so I protect it
 * from the auto-formatter that I use.
 */
// clang-format off
<% 
setup_pybind11(cfg) 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'].extend([ fp11.get_includes(), fp11.get_fwdpp_includes()])
#On OS X using clang, there is more work to do.  Using gcc on OS X
#gets rid of these requirements. The specifics sadly depend on how
#you initially built fwdpy11, and what is below assumes you used
#the provided setup.py + OS X + clang:
#cfg['compiler_args'].extend(['-stdlib=libc++','-mmacosx-version-min=10.7'])
#cfg['linker_args']=['-stdlib=libc++','-mmacosx-version-min=10.7']
#An alternative to the above is to add the first line to CPPFLAGS
#and the second to LDFLAGS when compiling a plugin on OS X using clang.
%>
// clang-format on

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>

    namespace py = pybind11;

struct snowdrift : public fwdpy11::SlocusPopGeneticValue
/* This is our stateful fitness object.
 * It records the model parameters and holds a
 * vector to track individual phenotypes.
 *
 * Here, we publicly inherit from fwdpy11::SlocusPopGeneticValue,
 * which is defined in the header included above.  It is
 * an abstract class in C++ terms, and is reflected
 * as a Python Abstract Base Class (ABC) called
 * fwdpy11.genetic_values.SlocusPopGeneticValue.
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

    snowdrift(double b1_, double b2_, double c1_, double c2_)
        : fwdpy11::SlocusPopGeneticValue{}, b1(b1_), b2(b2_), c1(c1_), c2(c2_),
          phenotypes()
    {
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop & /*pop*/) const
    // The call operator must return the genetic value of an individual
    {
        return phenotypes[diploid_index];
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
          const fwdpy11::SlocusPop & /*pop*/) const
    // This function may be used to model random effects...
    {
        //...but there are no random effects here.
        return 0.0;
    }

    inline void
    update(const fwdpy11::SlocusPop &pop)
    // A stateful fitness model needs updating.
    {
        phenotypes.resize(pop.N);
        for (std::size_t i = 0; i < pop.N; ++i)
            {
                // A diploid tracks its index via
                // fwdpy11::DiploidMetadata::label
                phenotypes[pop.diploid_metadata[i].label]
                    = fwdpp::additive_diploid(2.0)(pop.diploids[i],
                                                   pop.gametes, pop.mutations);
            }
    }
};

PYBIND11_MODULE(snowdrift, m)
{
    m.doc() = "Example of custom stateful fitness model.";

    // We need to import the Python version of our base class:
    pybind11::object imported_snowdrift_base_class_type
        = pybind11::module::import("fwdpy11.genetic_values")
              .attr("SlocusPopGeneticValue");

    // Create a Python class based on our new type
    py::class_<snowdrift, fwdpy11::SlocusPopGeneticValue>(m, "SlocusSnowdrift")
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"))
        .def_readwrite("phenotypes", &snowdrift::phenotypes);
}
