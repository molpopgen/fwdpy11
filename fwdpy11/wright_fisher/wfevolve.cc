#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <functional>
#include <cmath>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <iostream>

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
	std::cout << "here!\n";
    pop.mutations.reserve(std::ceil(
        std::log(2 * N)
        * (4. * double(N) * (mu_neutral + mu_selected)
           + 0.667 * (4. * double(N) * (mu_neutral + mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
			std::cout << generation << '\n';
            double wbar = KTfwd::sample_diploid(
                rng.get(), 
				pop.gametes, pop.diploids, pop.mutations,pop.mcounts, 
				pop.N, N, mu_neutral + mu_selected,
                KTfwd::extensions::bind_dmm(
                    mmodel, pop.mutations, pop.mut_lookup, rng.get(),
                    mu_neutral, mu_selected, pop.generation),
                recmap,
				fitness, 
				pop.neutral, pop.selected, selfing_rate);
            pop.N = N;
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
            recorder(pop, generation);
        }
}

namespace py = pybind11;

PYBIND11_PLUGIN(wfevolve)
{
    py::module m("wfevolve", "example extending");

    m.def("evolve_singlepop_regions_cpp", &evolve_singlepop_regions_cpp);

    return m.ptr();
}
