#include <exception>
#include <functional>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/multilocus.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

struct aggregate_additive_fitness
{
    inline double
    operator()(const py::array_t<double>& g) const noexcept
    {
        auto s = g.size();
        return std::max(0., std::accumulate(g.data(), g.data() + s, 0.0)
                                - (s - 1));
    }
};

struct aggregate_mult_fitness
{
    inline double
    operator()(const py::array_t<double>& g) const noexcept
    {
        auto s = g.size();
        return std::max(0., std::accumulate(g.data(), g.data() + s, 1.0,
                                            std::multiplies<double>()));
    }
};

struct aggregate_additive_trait
{
    inline double
    operator()(const py::array_t<double>& g) const noexcept
    {
        return std::accumulate(g.data(), g.data() + g.size(), 0.0);
    }
};

struct aggregate_mult_trait
{
    inline double
    operator()(const py::array_t<double>& g) const noexcept
    {
        auto s = g.size();
        return std::accumulate(
                   g.data(), g.data() + s, 1.0,
                   [](double prod, double v) { return prod * (1. + v); })
               - 1.0;
    }
};

#define AGGREGATOR(CPPNAME, PYNAME, DOCSTRING)                                \
    py::class_<CPPNAME>(m, PYNAME, DOCSTRING)                                 \
        .def(py::init<>())                                                    \
        .def("__call__", [](const CPPNAME& agg,                               \
                            const py::array_t<double>& a) { return agg(a); }) \
        .def("__getstate__",                                                  \
             [](const CPPNAME& aggregator) {                                  \
                 return py::make_tuple(std::string("CPPNAME"));               \
             })                                                               \
        .def("__setstate__", [](CPPNAME& aggregator, py::tuple t) {           \
            std::string n = t[0].cast<std::string>();                         \
            if (n == "CPPNAME")                                               \
                {                                                             \
                    new (&aggregator) CPPNAME();                              \
                }                                                             \
            else                                                              \
                {                                                             \
                    throw std::invalid_argument(                              \
                        "incorrect cppname encountered for aggregator");      \
                }                                                             \
        });

using ff_vec = decltype(fwdpy11::multilocus_genetic_value::fitness_functions);

PYBIND11_MODULE(multilocus, m)
{
    m.doc() = "Types and functions specific to multi-locus/region simulations";

    py::class_<fwdpy11::multilocus_genetic_value>(m, "MultiLocusGeneticValue")
        .def(py::init<const ff_vec&>())
        .def("__call__",
             [](const fwdpy11::multilocus_genetic_value& m,
                const fwdpy11::multilocus_diploid_t& dip,
                const fwdpy11::multilocus_t& pop) {
                 if (m.size() != pop.diploids[0].size())
                     {
                         throw std::invalid_argument("number of fitness "
                                                     "callbacks does not "
                                                     "equal number of loci");
                     }
                 return m(dip, pop.gametes, pop.mutations);
             })
        .def_readonly("fitness_functions",
                      &fwdpy11::multilocus_genetic_value::fitness_functions,
                      "Read-only access to stored fitenss functions")
        .def("__len__",
             [](const fwdpy11::multilocus_genetic_value& m) {
                 return m.size();
             })
        .def("__getstate__",
             [](const fwdpy11::multilocus_genetic_value& mw) {
                 /* This is some crazy magic:
                  * This object contains a vector of unique_ptr
                  * to single-locus functions.  We will clone each
                  * element as a shared_ptr wrapping the underling
                  * type.
                  * This works b/c pybind11 knows about the C++
                  * inheritance
                  * hierarchy.
                  */
                 py::list rv;
                 for (auto&& f : mw.fitness_functions)
                     {
                         rv.append(f->clone_shared());
                     }
                 return rv;
             })
        .def("__setstate__",
             /* The back-conversion is just as interesting.
              * We can simply type cast the list back to the
              * C++
              * vector type ff_vec, which is a typedef
              * defined above.
              */
             [](fwdpy11::multilocus_genetic_value& mw, py::list l) {
                 new (&mw) fwdpy11::multilocus_genetic_value(l.cast<ff_vec>());
             });

    AGGREGATOR(aggregate_additive_fitness, "AggAddFitness",
               "Map genetic values from a multi-locus diploid to fitness "
               "under an additive model.");
    AGGREGATOR(aggregate_additive_trait, "AggAddTrait",
               "Map genetic values from a multi-locus diploid to trait value "
               "under an additive model.");
    AGGREGATOR(aggregate_mult_fitness, "AggMultFitness",
               "Map genetic values from a multi-locus diploid to fitness "
               "under an multiplicative model.");
    AGGREGATOR(aggregate_mult_trait, "AggMultTrait",
               "Map genetic values from a multi-locus diploid to trait value "
               "under an multiplicative model.");

    py::class_<fwdpy11::interlocus_rec>(m, "InterlocusRecombination")
        .def(py::init<double, std::underlying_type<fwdpy11::interlocus_rec::
                                                       RECMODEL>::type>(),
             "This class parameterizes interlocus recombination in multilocus "
             "simulations.  You do not make instances of this class directly. "
             " Rather, you make calls to one of "
             ":func:`fwdpy11.multilocus.poisson_rec` or "
             ":func:`fwdpy11.multilocus.binomial_rec`.")
        .def("__getstate__",
             [](const fwdpy11::interlocus_rec& ir) {
                 return py::make_tuple(ir.param, ir.get_model());
             })
        .def("__setstate__",
             [](fwdpy11::interlocus_rec& ir, py::tuple t) {
                 double d = t[0].cast<double>();
                 auto m = t[1].cast<fwdpy11::interlocus_rec::mtype>();
                 new (&ir) fwdpy11::interlocus_rec(d, m);
             })
        .def("__repr__", [](const fwdpy11::interlocus_rec& ir) {
            std::string rv = "multilocus.InterlocusRecombination(";
            rv += std::to_string(ir.param);
            rv += ',';
            rv += std::to_string(ir.m);
            rv += ')';
            return rv;
        });

    m.def("poisson_rec",
          [](const std::vector<double>& means) {
              py::list rv;
              for (auto&& r : means)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::POISSON));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param means: A list of mean values.

    :rtype: list

    :return: A list of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.fwdpy11_types.GSLrng` as argument.
    )delim");

    m.def("poisson_rec",
          [](const double mean) {
              return fwdpy11::interlocus_rec(
                  mean, fwdpy11::interlocus_rec::RECMODEL::POISSON);
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param mean: The mean of a Poisson process.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.fwdpy11_types.GSLrng` as argument.
    )delim");

    m.def("binomial_rec",
          [](const std::vector<double>& probs) {
              py::list rv;
              for (auto&& r : probs)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param probs: A list of genetic distance in cM/100.

    :rtype: list

    :return: A list of of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.fwdpy11_types.GSLrng` as argument.
    )delim");

    m.def("binomial_rec",
          [](const double prob) {
              return fwdpy11::interlocus_rec(
                  prob, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL);
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param prob: The genetic distance in cM/100.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.fwdpy11_types.GSLrng` as argument.
    )delim");

    m.def("poisson_rec",
          [](const fwdpy11::GSLrng_t& rng, const std::vector<double>& means) {
              auto w = PyErr_WarnEx(PyExc_DeprecationWarning,
                                    "this overload of poisson_rec is "
                                    "deprecated.  Please use the version that "
                                    "only takes a list of means.",
                                    0);
              py::list rv;
              for (auto&& r : means)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::POISSON));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param means: A list of mean values.

    :rtype: list

    :return: A list of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. deprecated:: 0.1.3
    )delim");

    m.def("poisson_rec",
          [](const fwdpy11::GSLrng_t& rng, const double mean) {
              auto w = PyErr_WarnEx(PyExc_DeprecationWarning,
                                    "this overload of poisson_rec is "
                                    "deprecated.  Please use the version that "
                                    "only takes a mean argument.",
                                    0);
              return fwdpy11::interlocus_rec(
                  mean, fwdpy11::interlocus_rec::RECMODEL::POISSON);
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param mean: The mean of a Poisson process.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. deprecated:: 0.1.3
    )delim");

    m.def("binomial_rec",
          [](const fwdpy11::GSLrng_t& rng, const std::vector<double>& probs) {
              auto w = PyErr_WarnEx(PyExc_DeprecationWarning,
                                    "this overload of binomial_rec is "
                                    "deprecated.  Please use the version that "
                                    "only takes a list of probabilities.",
                                    0);
              py::list rv;
              for (auto&& r : probs)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param probs: A list of genetic distance in cM/100.

    :rtype: list

    :return: A list of of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. deprecated:: 0.1.3
    )delim");

    m.def("binomial_rec",
          [](const fwdpy11::GSLrng_t& rng, const double prob) {
              auto w = PyErr_WarnEx(PyExc_DeprecationWarning,
                                    "this overload of binomial_rec is "
                                    "deprecated.  Please use the version that "
                                    "only takes a probability argument.",
                                    0);
              return fwdpy11::interlocus_rec(
                  prob, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL);
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param prob: The genetic distance in cM/100.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. deprecated:: 0.1.3
    )delim");
}
