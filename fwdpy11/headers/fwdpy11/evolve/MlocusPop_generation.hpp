//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef FWDPY11_EVOLVE_MLOCUSPOP_GENERATION_HPP__
#define FWDPY11_EVOLVE_MLOCUSPOP_GENERATION_HPP__

#include <tuple>
#include <vector>
#include <stdexcept>
#include <functional>
#include <type_traits>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/genetic_values/MlocusPopGeneticValue.hpp>
#include <gsl/gsl_randist.h>

namespace fwdpy11
{
    template <typename poptype, typename pick1_function,
              typename pick2_function, typename update_function,
              typename mutation_model, typename recombination_model>
    void
    evolve_generation(
        const GSLrng_t& rng, poptype& pop, const fwdpp::uint_t N_next,
        const std::vector<double>& mutrates, const mutation_model& mmodel,
        const recombination_model& recmodel,
        const std::vector<std::function<unsigned(void)>>& interlocus_rec,
        const pick1_function& pick1, const pick2_function& pick2,
        const update_function& update)
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   fwdpp::sugar::MULTILOC_TAG>::value,
                      "Population type must be a multi-locus type.");

        auto gamete_recycling_bin
            = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
        auto mutation_recycling_bin
            = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);

        // Efficiency hit.  Unavoidable
        // in use case of a sampler looking
        // at the gametes themselves (even tho
        // gamete.n has little bearing on anything
        // beyond recycling).  Can revisit later
        for (auto&& g : pop.gametes)
            g.n = 0;

        decltype(pop.diploids) offspring(N_next);
        decltype(pop.diploid_metadata) offspring_metadata(N_next);
        // Generate the offspring
        std::size_t label = 0;
        for (auto& dip : offspring)
            {
                auto p1 = pick1();
                auto p2 = pick2(p1);

                dip = fwdpp::fwdpp_internal::multilocus_rec_mut(
                    rng.get(), pop.diploids[p1], pop.diploids[p2],
                    mutation_recycling_bin, gamete_recycling_bin, recmodel,
                    interlocus_rec,
                    ((gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0),
                    ((gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0), pop.gametes,
                    pop.mutations, pop.neutral, pop.selected, mutrates.data(),
                    mmodel);

#ifndef NDEBUG
                for (const auto& locus : dip)
                    {
                        if (pop.gametes[locus.first].n == 0
                            || pop.gametes[locus.second].n == 0)
                            {
                                throw std::runtime_error(
                                    "DEBUG: locus has gamete with frequency zero");
                            }
                    }
#endif
                offspring_metadata[label].label = label;
                update(offspring_metadata[label++], p1, p2,
                       pop.diploid_metadata);
            }

        fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                               pop.mcounts);
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
    }
} // namespace fwdpy11

#endif
