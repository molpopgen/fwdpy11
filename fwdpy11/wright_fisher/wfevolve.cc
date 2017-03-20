#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <functional>
#include <cmath>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include "../src/types.hpp"

// Evolve the population for some amount of time with mutation and recombination
void evolve(
    fwdpy::singlepop_t& pop, const fwdpy::GSLrng_t & rng,  const unsigned& N,
    const unsigned& generations, const double& mu, const double& recrate,
    std::function<void(const fwdpy::singlepop_t&, const unsigned)> recorder) {
    pop.mutations.reserve(std::ceil(std::log(2 * N) * (4. * double(N) * (mu)) +
                                    0.667 * (4. * double(N) * (mu))));
    std::function<double(void)> recmap =
        std::bind(gsl_rng_uniform, rng.get());  // uniform crossover map
    for (unsigned generation = 0; generation < generations; ++generation) {
        double wbar = KTfwd::sample_diploid(
            rng.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts, N,
            mu, std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, rng.get(),
                          std::ref(pop.mut_lookup), generation,mu, 0.,
                          [&rng]() { return gsl_rng_uniform(rng.get()); },
                          []() { return 0.; }, []() { return 0.; }),
            std::bind(KTfwd::poisson_xover(), rng.get(), recrate, 0., 1.,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                      std::placeholders::_2, std::placeholders::_3, 2.),
            pop.neutral, pop.selected);
        KTfwd::update_mutations(pop.mutations, pop.fixations,
                                pop.fixation_times, pop.mut_lookup, pop.mcounts,
                                generation, 2 * pop.N);
        recorder(pop, generation);
    }
}

namespace py = pybind11;

PYBIND11_PLUGIN(wfevolve) {
    py::module m("wfevolve", "example extending");

    m.def("evolve", &evolve);

    return m.ptr();
}
