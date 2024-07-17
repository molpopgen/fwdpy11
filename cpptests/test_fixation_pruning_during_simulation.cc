#include <algorithm>
#include <boost/test/tools/fpc_op.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/unit_test_suite.hpp>
#include <core/evolve_discrete_demes/evolvets.hpp>
#include <core/demes/forward_graph.hpp>
#include <cstdint>
#include <cwchar>
#include <exception>
#include <fwdpy11/evolvets/SampleRecorder.hpp>
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/genetic_values/fwdpp_wrappers/fwdpp_genetic_value.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/genetic_values/DiploidMultiplicative.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/regions/ExpS.hpp>
#include <memory>
#include <stdexcept>

#include "forward_demes_graph_fixtures.hpp"
#include "fwdpp/ts/simplification/simplification.hpp"
#include "fwdpy11/discrete_demography/exceptions.hpp"
#include "fwdpy11/genetic_values/DiploidMultiplicative.hpp"
#include "fwdpy11/mutation_dominance/MutationDominance.hpp"
#include "fwdpy11/regions/Sregion.hpp"

fwdpy11::MutationRegions
strong_positive_selection()
{
    auto nweights = std::vector<double>();
    auto sweights = std::vector<double>(1.);
    std::vector<std::unique_ptr<fwdpy11::Sregion>> nregions, sregions;

    sregions.emplace_back(fwdpy11::ExpS(fwdpy11::Region(0., 10., 1.0, true, 0), 2.0,
                                        0.025, fwdpy11::fixed_dominance(1.))
                              .clone());
    return fwdpy11::MutationRegions::create(0., nweights, sweights, nregions, sregions);
}

struct common_setup
{
    fwdpy11::GSLrng_t rng;
    fwdpy11::DiploidPopulation pop;
    fwdpy11::MutationRegions mregions;
    fwdpy11::RecombinationRegions recregions;
    fwdpy11::DiploidMultiplicative multiplicative;
    fwdpy11::dgvalue_pointer_vector_ gvalue_ptrs;
    fwdpy11::SampleRecorder recorder;
    fwdpy11_core::ForwardDemesGraph forward_demes_graph;
    std::function<void(const fwdpy11::DiploidPopulation &)> post_simplification_recorder;
    std::function<bool(const fwdpy11::DiploidPopulation &, const bool)>
        stopping_criterion;
    evolve_with_tree_sequences_options options;

    common_setup()
        : rng{42}, pop{100, 10.0}, mregions{strong_positive_selection()},
          recregions{0., {}}, multiplicative{1,
                                             2.0,
                                             fwdpy11::final_multiplicative_fitness(),
                                             [](const auto) { return false; },
                                             nullptr,
                                             nullptr},
          gvalue_ptrs{multiplicative}, recorder{},
          forward_demes_graph(SingleDemeModel().yaml, 10000),
          post_simplification_recorder{[](const fwdpy11::DiploidPopulation &) {}},
          stopping_criterion{[](const fwdpy11::DiploidPopulation &, const bool) -> bool {
              return false;
          }},
          options{}
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_remove_fixations_multiplicative)

BOOST_FIXTURE_TEST_CASE(test_remove_fixations, common_setup)
{
    const auto sample_recorder_callback = [](const fwdpy11::DiploidPopulation &pop,
                                             fwdpy11::SampleRecorder &) {
        for (const auto diploid : pop.diploids)
            {
                for (const auto m : pop.haploid_genomes[diploid.first].smutations)
                    {
                        if (pop.mcounts[m] == 2 * pop.N)
                            {
                                for (std::size_t i = 0; i < pop.fixations.size(); ++i)
                                    {
                                        if (pop.fixations[i].g == pop.mutations[m].g
                                            && pop.fixations[i].pos
                                                   == pop.mutations[m].pos
                                            && pop.fixations[i].s == pop.mutations[m].s)
                                            {
                                                if (pop.fixation_times[i]
                                                    != pop.generation)
                                                    {
                                                        throw std::runtime_error(
                                                            "fixation found in genome "
                                                            "that did not occur this "
                                                            "generation");
                                                    }
                                            }
                                    }
                            }
                    }
            }
    };
    BOOST_REQUIRE_EQUAL(mregions.weights.size(), mregions.regions.size());
    BOOST_REQUIRE_EQUAL(mregions.weights.size(), 1);
    options.track_mutation_counts_during_sim = true;
    BOOST_REQUIRE_EQUAL(options.preserve_selected_fixations, false);
    BOOST_REQUIRE_EQUAL(options.track_mutation_counts_during_sim, true);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 1000, 0.,
                               0.02, mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 1000);
    BOOST_REQUIRE_EQUAL(pop.fixations.size(), pop.fixation_times.size());
    BOOST_REQUIRE(!pop.fixations.empty());
}

BOOST_AUTO_TEST_SUITE_END()
