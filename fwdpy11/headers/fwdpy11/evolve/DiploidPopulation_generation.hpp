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
#ifndef FWDPY11_EVOLVE_POP_GENERATION_HPP__
#define FWDPY11_EVOLVE_POP_GENERATION_HPP__

#include <tuple>
#include <type_traits>
#include <stdexcept>
#include <fwdpp/internal/haploid_genome_cleaner.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>
#include <gsl/gsl_randist.h>

namespace fwdpy11
{
    template <typename poptype, typename pick1_function,
              typename pick2_function, typename update_function,
              typename mutation_model, typename recombination_model>
    void
    evolve_generation(const GSLrng_t& rng, poptype& pop,
                      const fwdpp::uint_t N_next, const double mu,
                      const mutation_model& mmodel,
                      const recombination_model& recmodel,
                      const pick1_function& pick1, const pick2_function& pick2,
                      const update_function& update)
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   fwdpp::poptypes::DIPLOID_TAG>::value,
                      "Population type must be a diploid type.");

        auto gamete_recycling_bin
            = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);
        auto mutation_recycling_bin = fwdpp::make_mut_queue(pop.mcounts);

        // Efficiency hit.  Unavoidable
        // in use case of a sampler looking
        // at the haploid_genomes themselves (even tho
        // gamete.n has little bearing on anything
        // beyond recycling).  Can revisit later
        for (auto&& g : pop.haploid_genomes)
            g.n = 0;

        decltype(pop.diploids) offspring(N_next);
        decltype(pop.diploid_metadata) offspring_metadata(N_next);
        // Generate the offspring
        std::size_t label = 0;
        for (auto& dip : offspring)
            {
                auto p1 = pick1();
                auto p2 = pick2(p1);

                auto p1g1 = pop.diploids[p1].first;
                auto p1g2 = pop.diploids[p1].second;
                auto p2g1 = pop.diploids[p2].first;
                auto p2g2 = pop.diploids[p2].second;

                // Mendel
                if (gsl_rng_uniform(rng.get()) < 0.5)
                    std::swap(p1g1, p1g2);
                if (gsl_rng_uniform(rng.get()) < 0.5)
                    std::swap(p2g1, p2g2);

                // Update to fwdpp 0.5.7 API
                // in fwdpy11 0.1.4
                fwdpp::mutate_recombine_update(
                    rng.get(), pop.haploid_genomes, pop.mutations,
                    std::make_tuple(p1g1, p1g2, p2g1, p2g2), recmodel, mmodel,
                    mu, gamete_recycling_bin, mutation_recycling_bin, dip,
                    pop.neutral, pop.selected);

#ifndef NDEBUG
                if (pop.haploid_genomes[dip.first].n == 0
                    || pop.haploid_genomes[dip.second].n == 0)
                    {
                        throw std::runtime_error(
                            "DEBUG: diploid has gamete with frequency zero");
                    }
#endif
                offspring_metadata[label].label = label;
                update(offspring_metadata[label++], p1, p2,
                       pop.diploid_metadata);
            }

        fwdpp::fwdpp_internal::process_haploid_genomes(
            pop.haploid_genomes, pop.mutations, pop.mcounts);
        // This is constant-time
        pop.diploids.swap(offspring);
        pop.diploid_metadata.swap(offspring_metadata);
    }
} // namespace fwdpy11

#endif
