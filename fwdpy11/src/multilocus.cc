#include <exception>
#include <functional>
#include <algorithm>
#include <limits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpy11/rng.hpp>
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
        .def("__call__",                                                      \
             [](const CPPNAME& agg, const py::array_t<double>& a) {           \
                 return agg(a);                                               \
             });

using ff_vec = std::vector<std::shared_ptr<fwdpy11::single_locus_fitness>>;

PYBIND11_PLUGIN(multilocus)
{
    py::module m(
        "multilocus",
        "Types and functions specific to multi-locus/region simulations");

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
        .def("__len__", [](const fwdpy11::multilocus_genetic_value& m) {
            return m.size();
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

    m.def("poisson_rec",
          [](const fwdpy11::GSLrng_t& rng, const std::vector<double>& rates) {
              return KTfwd::make_poisson_interlocus_rec(rng.get(), rates.data(),rates.size());
          });

    m.def("poisson_rec",
          [](const fwdpy11::GSLrng_t& rng, const double rate) -> std::function<unsigned()> {
              return std::bind(gsl_ran_poisson,rng.get(),rate); 
          });

    m.def("binomial_rec",
          [](const fwdpy11::GSLrng_t& rng, const std::vector<double>& rates) {
              return KTfwd::make_binomial_interlocus_rec(rng.get(), rates.data(),rates.size());
          });

    m.def("binomial_rec",
          [](const fwdpy11::GSLrng_t& rng, const double prob) -> std::function<unsigned()> {
              return std::bind(gsl_ran_binomial,rng.get(),prob,1); 
          });

    return m.ptr();
}
