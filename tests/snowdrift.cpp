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

struct snowdrift_diploid
/* This is a function object implementing the snowdrift
 * fitness calculation
 */
{
    using result_type = double;
    inline result_type
    operator()(const fwdpy11::DiploidGenotype &dip,
               const fwdpy11::Population::gcont_t &gametes,
               const fwdpy11::Population::mcont_t &mutations,
               const fwdpy11::DiploidMetadata&metadata,
               const std::vector<double> &phenotypes, const double b1,
               const double b2, const double c1, const double c2) const
    /* The first 3 arguments will be passed in from fwdpp.
     * The phenotypes are the stateful part, and will get
     * passed in by reference.  The remaining params
     * define the payoff function and are also supplied externally.
     */
    {
        auto N = phenotypes.size();
        // A diploid tracks its index via
        // fwdpy11::Diploid::label, which
        // is std::size_t.
        auto i = metadata.label;
        double zself = phenotypes[i];
        // The code here would not be found
        // in production code.  Here, we are testing
        // that all the data have been updated correctly
        // in fwdpy11's machinery for evolving populations.
        // This test is only here as this is part of unit
        // testing suite.
        auto g = fwdpp::additive_diploid(2.0)(dip, gametes, mutations);
        if (zself != g)
            {
                throw std::runtime_error("snowdrift not working: "
                                         + std::to_string(zself) + ' '
                                         + std::to_string(g));
            }
        // End unit-testing only code
        result_type fitness = 0;
        for (std::size_t j = 0; j < N; ++j)
            {
                if (i != j)
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
};

struct snowdrift : public fwdpy11::SlocusPopGeneticValue
/* This is our stateful fitness object.
 * It records the model parameters and holds a
 * vector to track individual phenotypes.
 *
 * The C++ side of an SlocusFitness object must publicly
 * inherit from fwdpy11::single_locus_fitness.
 *
 * The phenotypes get updated each generation during
 * the simulation.
 *
 * The phenotypes will be the simple additive model,
 * calculated using fwdpp's machinery.
 */
{
    double b1, b2, c1, c2;
    std::vector<double> phenotypes;

    snowdrift(double b1_, double b2_, double c1_, double c2_)
        : fwdpy11::SlocusPopGeneticValue{}, b1(b1_), b2(b2_), c1(c1_), c2(c2_),
          phenotypes()
    {
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop &pop) const
    {
        return fwdpp::additive_diploid(2.0)(pop.diploids[diploid_index],
                                            pop.gametes, pop.mutations);
    }

    inline double
    genetic_value_to_fitness(const fwdpy11::DiploidMetadata &metadata) const
    {
        double fitness = 0.0;
        double zself = metadata.g;
        auto N = phenotypes.size();
        for (std::size_t j = 0; j < N; ++j)
            {
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
          const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
          const std::size_t /*parent1*/, const std::size_t /*parent2*/,
          const fwdpy11::SlocusPop & /*pop*/) const
    {
        return 0.0;
    }

    inline void
    update(const fwdpy11::SlocusPop &pop)
    /* A stateful fitness model needs updating.
     * The base class defines this virtual function
     * to do nothing (for non-stateful models).
     * Here, we redefine it as needed.
     */
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

    //pybind11::object FWDPY11_SINGLE_LOCUS_FITNESS_BASE_IMPORT__               \
    //    = (pybind11::object)pybind11::module::import("fwdpy11.fitness")       \
    //          .attr("SlocusFitness");
    try
        {
            pybind11::object imported_base
                = pybind11::module::import("fwdpy11.genetic_values")
                      .attr("fwdpy11.genetic_values.SlocusPopGeneticValue");
        }
    catch (...)
        {
        }
    // Create a Python class based on our new type
    py::class_<snowdrift, fwdpy11::SlocusPopGeneticValue>(m, "SlocusSnowdrift")
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"))
        .def_readwrite("phenotypes", &snowdrift::phenotypes);
    // It is useful to make these type compatible with Python's
    // pickling protocol:
    //.def(py::pickle(
    //    [](const snowdrift &s) {
    //        return py::make_tuple(s.b1, s.b2, s.c1, s.c2, s.phenotypes);
    //    },
    //    [](py::tuple t) {
    //        auto rv = std::make_shared<snowdrift>(
    //            t[0].cast<double>(), t[1].cast<double>(),
    //            t[2].cast<double>(), t[3].cast<double>());
    //        rv->phenotypes = t[4].cast<std::vector<double>>();
    //        return rv;
    //    }));
}
