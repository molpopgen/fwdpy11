#include <fwdpy11/types.hpp>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <functional>
#include <cmath>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/fitness/fitness.hpp>

// Evolve the population for some amount of time with mutation and
// recombination
void
evolve_singlepop_regions_cpp(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop, const unsigned& N,
    const unsigned generations, const double mu_neutral,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    const fwdpy11::singlepop_fitness& fitness,
    fwdpy11::singlepop_temporal_sampler recorder, const double selfing_rate)
{
    const auto fitness_callback = fitness.callback();
    pop.mutations.reserve(std::ceil(
        std::log(2 * N)
        * (4. * double(N) * (mu_neutral + mu_selected)
           + 0.667 * (4. * double(N) * (mu_neutral + mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            double wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, N, mu_neutral + mu_selected,
                KTfwd::extensions::bind_dmm(
                    mmodel, pop.mutations, pop.mut_lookup, rng.get(),
                    mu_neutral, mu_selected, pop.generation),
                recmap, fitness_callback,
                pop.neutral, pop.selected, selfing_rate);
            pop.N = N;
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
            recorder(pop, generation);
        }
}

// Evolve the population for some amount of time with mutation and
// recombination
void
evolve(fwdpy11::singlepop_t& pop, const fwdpy11::GSLrng_t& rng,
       const unsigned& N, const unsigned& generations, const double& mu,
       const double& recrate)
{
    pop.mutations.reserve(std::ceil(std::log(2 * N) * (4. * double(N) * (mu))
                                    + 0.667 * (4. * double(N) * (mu))));
    for (unsigned generation = 0; generation < generations; ++generation)
        {
            double wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, N, mu + 0.005,
                std::bind(
                    KTfwd::infsites(), std::placeholders::_1,
                    std::placeholders::_2, rng.get(), std::ref(pop.mut_lookup),
                    &generation, mu, 0.005,
                    [&rng]() { return gsl_rng_uniform(rng.get()); },
                    [&rng]() { return gsl_ran_exponential(rng.get(), -0.1); },
                    []() { return 1.; }),
                std::bind(KTfwd::poisson_xover(), rng.get(), recrate, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3, 2.),
                pop.neutral, pop.selected);
            pop.N = N;
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
        }
}

namespace py = pybind11;

PYBIND11_PLUGIN(wfevolve)
{
    py::module m("wfevolve", "example extending");

    m.def("evolve_singlepop_regions_cpp", &evolve_singlepop_regions_cpp);
    m.def("evolve", &evolve);
    return m.ptr();
}
