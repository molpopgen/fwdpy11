#ifndef FWDPY11_EVOLVE_MLOCUSPOP_HPP__
#define FWDPY11_EVOLVE_MLOCUSPOP_HPP__

#include <tuple>
#include <type_traits>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/recombination.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/samplers.hpp>
#include <gsl/gsl_randist.h>

namespace fwdpy11
{
    template <typename poptype, typename pick1_function,
              typename pick2_function, typename update_function,
              typename mutation_model, typename recombination_model,
              typename mutation_removal_policy>
    void
    evolve_generation(
        const GSLrng_t& rng, poptype& pop, const KTfwd::uint_t N_next,
        const std::vector<double>& mu, const mutation_model& mmodel,
        const recombination_model& recmodel,
        const std::vector<std::function<unsigned(void)>>& interlocus_rec,
        const pick1_function& pick1, const pick2_function& pick2,
        const update_function& update, const mutation_removal_policy& mrp)
    {
        static_assert(
            std::is_same<typename poptype::popmodel_t,
                         KTfwd::sugar::MULTILOCPOP_TAG>::value,
            "Population type must be a multi-locus, single-deme type.");

        auto gamete_recycling_bin
            = KTfwd::fwdpp_internal::make_gamete_queue(pop.gametes);
        auto mutation_recycling_bin
            = KTfwd::fwdpp_internal::make_mut_queue(pop.mcounts);

        for (auto&& g : pop.gametes)
            g.n = 0;

        decltype(pop.diploids) offspring(N_next);

        // Generate the offspring
		std::size_t label = 0;
        for (auto& dip : offspring)
            {
                auto p1 = pick1(rng, pop);
                auto p2 = pick2(rng, pop, p1);

                dip = KTfwd::fwdpp_internal::multilocus_rec_mut(
                    rng.get(), pop.diploids[p1], pop.diploids[p2],
                    mutation_recycling_bin, gamete_recycling_bin, recmodel,
                    interlocus_rec,
                    ((gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0),
                    ((gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0), pop.gametes,
                    pop.mutations, pop.neutral, pop.selected, mu.data(),
                    mmodel, KTfwd::emplace_back());
				dip[0].label=label++;
                update(rng, dip, pop, p1, p2);
            }

        KTfwd::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                               pop.mcounts);
        KTfwd::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                              pop.mcounts, 2 * N_next, mrp,
                                              std::true_type());
        // This is constant-time
        pop.diploids.swap(offspring);
    }
}

#endif
