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
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

    namespace py = pybind11;

struct snowdrift_diploid
/* This is a function object implementing the snowdrift
 * fitness calculation
 */
{
    using result_type = double;
    inline result_type
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations,
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
        // fwdpy11::diploid_t::label, which
        // is std::size_t.
        auto i = dip.label;
        double zself = phenotypes[i];
        // The code here would not be found
        // in production code.  Here, we are testing
        // that all the data have been updated correctly
        // in fwdpy11's machinery for evolving populations.
        // This test is only here as this is part of unit
        // testing suite.
        auto g = KTfwd::additive_diploid()(dip, gametes, mutations, 2.0);
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

struct snowdrift : public fwdpy11::single_locus_fitness
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
        : b1(b1_), b2(b2_), c1(c1_), c2(c2_), phenotypes(std::vector<double>())
    {
    }

    inline fwdpy11::single_locus_fitness_fxn
    callback() const final
    // A stateful fitness model must return a bound callable.
    {
        /* Please remember that std::bind will pass a copy unless you
         * explicity wrap with a reference wrapper.  Here,
         * we use std::cref to make sure that phenotypes are
         * passed as const reference.
         */
        return std::bind(snowdrift_diploid(), std::placeholders::_1,
                         std::placeholders::_2, std::placeholders::_3,
                         std::cref(phenotypes), b1, b2, c1, c2);
    }
    void
    update(const fwdpy11::singlepop_t &pop) final
    /* A stateful fitness model needs updating.
     * The base class defines this virtual function
     * to do nothing (for non-stateful models).
     * Here, we redefine it as needed.
     */
    {
        phenotypes.resize(pop.N);
        for (auto &&dip : pop.diploids)
            {
                // A diploid tracks its index via
                // fwdpy11::diploid_t::label
                phenotypes[dip.label] = KTfwd::additive_diploid()(
                    dip, pop.gametes, pop.mutations, 2.0);
            }
    };

    // A custom fitness function requires that
    // several pure virtual functions be defined.
    // These macros do it for you:
    SINGLE_LOCUS_FITNESS_CLONE_SHARED(snowdrift);
    SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(snowdrift);
    // This one gets pass the C++ name of the
    // callback.  We use C++11's typeid to
    // do that for us:
    SINGLE_LOCUS_FITNESS_CALLBACK_NAME(typeid(snowdrift_diploid).name());
};

PYBIND11_MODULE(snowdrift, m)
{
    m.doc() = "Example of custom stateful fitness model.";

    // fwdpy11 provides a macro
    // to make sure that the Python wrapper
    // to fwdpy11::singlepop_fitness is visible
    // to this module.
    FWDPY11_SINGLE_LOCUS_FITNESS()

    // Create a Python class based on our new type
    py::class_<snowdrift, std::shared_ptr<snowdrift>,
               fwdpy11::single_locus_fitness>(m, "SlocusSnowdrift")
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"))
        .def_readwrite("phenotypes", &snowdrift::phenotypes)
        // It is useful to make these type compatible with Python's
        // pickling protocol:
        .def("__getstate__",
             [](const snowdrift &s) {
                 return py::make_tuple(s.b1, s.b2, s.c1, s.c2, s.phenotypes);
             })
        .def("__setstate__", [](snowdrift &s, py::tuple t) {
            new (&s) snowdrift(t[0].cast<double>(), t[1].cast<double>(),
                               t[2].cast<double>(), t[3].cast<double>());
            s.phenotypes = t[4].cast<std::vector<double>>();
        });
}
