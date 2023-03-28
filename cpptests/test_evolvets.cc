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
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <stdexcept>

#include "forward_demes_graph_fixtures.hpp"
#include "fwdpp/ts/simplification/simplification.hpp"
#include "fwdpy11/discrete_demography/exceptions.hpp"

/* Tests needed
 *
 * - [X] simlen > model time in graph
 * - [X] initial population config not compatible
 *       with state of parental demes for single deme models.
 * - [X] initial population config not compatible
 *       with state of parental demes for multi deme models.
 * - [X] simlen < model time in graph,
 *       but we keep simulating until we are done.
 * - [X] One deme, multiple epochs, test size history is 
 *      correct
 */

namespace
{
    // NOTE: all of this is a copy/paste
    // from code internal to the pybind11 module.

    struct single_deme_additive_het
    {
        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += m.s * m.h;
        }
    };

    struct multi_deme_additive_het
    {
        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += m.esizes[deme] * m.heffects[deme];
        }
    };

    struct single_deme_additive_hom
    {
        double scaling;
        single_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.s;
        }
    };

    struct multi_deme_additive_hom
    {
        double scaling;
        multi_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.esizes[deme];
        }
    };

    struct final_additive_trait
    {
        inline double
        operator()(double d) const
        {
            return d;
        }
    };

    struct final_additive_fitness
    {
        inline double
        operator()(double d) const
        {
            return std::max(0.0, 1. + d);
        }
    };

    using DiploidAdditive = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
        single_deme_additive_het, single_deme_additive_hom, multi_deme_additive_het,
        multi_deme_additive_hom, 0>;
}

template <typename T> struct common_setup
{
    fwdpy11::GSLrng_t rng;
    fwdpy11::DiploidPopulation pop;
    fwdpy11::MutationRegions mregions;
    fwdpy11::RecombinationRegions recregions;
    DiploidAdditive additive;
    fwdpy11::dgvalue_pointer_vector_ gvalue_ptrs;
    fwdpy11::SampleRecorder recorder;
    std::function<void(const fwdpy11::DiploidPopulation&)> post_simplification_recorder;
    std::function<void(const fwdpy11::DiploidPopulation&, fwdpy11::SampleRecorder&)>
        sample_recorder_callback;
    std::function<bool(const fwdpy11::DiploidPopulation&, const bool)>
        stopping_criterion;
    evolve_with_tree_sequences_options options;

    common_setup()
        : rng{42}, pop{100, 10.0}, mregions{{}, {}}, recregions{0., {}},
          additive{1, 2.0, final_additive_fitness(), nullptr, nullptr},
          gvalue_ptrs{additive}, recorder{},
          post_simplification_recorder{[](const fwdpy11::DiploidPopulation&) {}},
          sample_recorder_callback{
              [](const fwdpy11::DiploidPopulation&, fwdpy11::SampleRecorder&) {}},
          stopping_criterion{[](const fwdpy11::DiploidPopulation&, const bool) -> bool {
              return false;
          }},
          options{}
    {
    }

    T
    build_model() const
    {
        return T();
    }
};

template <typename T>
struct common_setup_with_sample_history_recording : public common_setup<T>
{
    std::unordered_map<std::int32_t,
                       std::vector<std::pair<std::uint32_t, std::uint32_t>>>
        size_history;
    common_setup_with_sample_history_recording() : common_setup<T>(), size_history{}
    {
        this->sample_recorder_callback = [this](const fwdpy11::DiploidPopulation& pop,
                                                fwdpy11::SampleRecorder&) {
            std::unordered_map<std::int32_t, std::uint32_t> deme_sizes;
            for (auto& md : pop.diploid_metadata)
                {
                    if (deme_sizes.find(md.deme) != std::end(deme_sizes))
                        {
                            deme_sizes[md.deme]++;
                        }
                    else
                        {
                            deme_sizes.emplace(md.deme, 1);
                        }
                }
            for (const auto& x : deme_sizes)
                {
                    if (this->size_history.find(x.first) == std::end(this->size_history))
                        {
                            std::vector<std::pair<std::uint32_t, std::uint32_t>> h{
                                std::make_pair(pop.generation, x.second)};
                            this->size_history.emplace(x.first, std::move(h));
                        }
                    else
                        {
                            this->size_history[x.first].push_back(
                                std::make_pair(pop.generation, x.second));
                        }
                }
        };
    }
};

struct ancestry_proportions
{
    std::uint32_t generation;
    std::int32_t offspring_deme;
    std::int32_t parent_deme;
    std::uint32_t num_parents;
};

// NOTE: efficiency is out the window here.
template <typename T> struct common_setup_with_ancestry_tracking : public common_setup<T>
{
    std::vector<ancestry_proportions> ancestry;
    // parent id -> parent deme
    std::unordered_map<std::size_t, std::int32_t> parents;
    common_setup_with_ancestry_tracking() : common_setup<T>(), ancestry(), parents()
    {
        this->sample_recorder_callback = [this](const fwdpy11::DiploidPopulation& pop,
                                                fwdpy11::SampleRecorder&) {
            std::vector<ancestry_proportions> current_ancestry;
            for (auto& md : pop.diploid_metadata)
                {
                    for (auto p : md.parents)
                        {
                            if (parents.find(p) == std::end(parents))
                                {
                                    throw std::runtime_error("parent not found");
                                }
                        }
                    if (parents[md.parents[0]] != parents[md.parents[1]])
                        {
                            throw std::runtime_error("parents from different demes");
                        }
                    auto parent_deme = parents[md.parents[0]];
                    auto iter = std::find_if(std::begin(current_ancestry),
                                             std::end(current_ancestry),
                                             [md, parent_deme](const auto& x) {
                                                 return x.offspring_deme == md.deme
                                                        && x.parent_deme == parent_deme;
                                             });
                    if (iter == std::end(current_ancestry))
                        {
                            current_ancestry.push_back(ancestry_proportions{
                                pop.generation, md.deme, parent_deme, 1});
                        }
                    else
                        {
                            iter->num_parents++;
                        }
                }
            for (auto& a : current_ancestry)
                {
                    ancestry.push_back(a);
                }
            this->update_parents(pop);
        };
    }

    void
    update_parents(const fwdpy11::DiploidPopulation& pop)
    {
        parents.clear();
        for (auto& md : pop.diploid_metadata)
            {
                parents.emplace(md.label, md.deme);
            }
    }
};

bool
validate_ancestry(const std::vector<ancestry_proportions>& ancestry,
                  std::uint32_t generation, std::int32_t offspring_deme,
                  std::int32_t parent_deme, std::function<bool(std::uint32_t)> condition)
{
    for (auto& a : ancestry)
        {
            if (a.generation == generation && a.offspring_deme == offspring_deme
                && a.parent_deme == parent_deme)
                {
                    BOOST_REQUIRE(condition(a.num_parents));
                    return true;
                }
        }
    return false;
}

BOOST_AUTO_TEST_SUITE(test_evolvets)

BOOST_FIXTURE_TEST_CASE(test_basic_api_coherence, common_setup<SingleDemeModel>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // TODO: if we put long run times in here, we get exceptions
    // from the ForwardDemesGraph back end.
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 10, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_FIXTURE_TEST_CASE(test_basic_api_coherence_two_deme_perpetual_island_model,
                        common_setup<TwoDemePerpetualIslandModel>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    // over-write the fixture so that the initial pop is okay
    pop = fwdpy11::DiploidPopulation(
        forward_demes_graph.parental_deme_sizes_at_time_zero(), 10.0);

    // TODO: if we put long run times in here, we get exceptions
    // from the ForwardDemesGraph back end.
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 10, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
    std::vector<unsigned> ndemes{0, 0};
    for (const auto& dip : pop.diploid_metadata)
        {
            ndemes[dip.deme]++;
        }
    BOOST_REQUIRE_EQUAL(ndemes[0], 100);
    BOOST_REQUIRE_EQUAL(ndemes[1], 100);
}

BOOST_FIXTURE_TEST_CASE(test_simlen_longer_than_model_length,
                        common_setup<SingleDemeModel>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // Here, simlen = 100, but the graph only has 10 generations of births in it.
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 100, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_FIXTURE_TEST_CASE(test_invalid_deme_metadata, common_setup<SingleDemeModel>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // Make some of the individuals in deme 1 but all must be in deme 0
    // b/c that is the model
    for (std::size_t i = 0; i < pop.diploid_metadata.size(); ++i)
        {
            if (i % 2 == 0)
                {
                    pop.diploid_metadata[i].deme = 1;
                }
        }

    BOOST_CHECK_THROW(
        {
            evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 10,
                                       0., 0., mregions, recregions, gvalue_ptrs,
                                       sample_recorder_callback, stopping_criterion,
                                       post_simplification_recorder, options);
        },
        fwdpy11::discrete_demography::DemographyError);
    // pop hasn't evolved!
    BOOST_REQUIRE_EQUAL(pop.generation, 0);
}

BOOST_FIXTURE_TEST_CASE(test_initial_pop_size_invalid, common_setup<SingleDemeModel>)
{
    // The demes model specifies N = 100.
    // Here, we will start with a different N.
    // This is a hard error!
    // We test values ~100 because TDD found cases of memory errors
    // when N >> the correct N but things can pass when N =~ the
    // correct N. Gotta love UB :).
    std::vector<std::uint32_t> initial_popsizes{50, 99, 101, 200};
    for (auto initial_n : initial_popsizes)
        {
            // reset the fixture
            pop = fwdpy11::DiploidPopulation(initial_n, 10.0);
            auto model = build_model();
            fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

            BOOST_CHECK_THROW(
                {
                    evolve_with_tree_sequences(
                        rng, pop, recorder, 10, forward_demes_graph, 10, 0., 0.,
                        mregions, recregions, gvalue_ptrs, sample_recorder_callback,
                        stopping_criterion, post_simplification_recorder, options);
                },
                // TODO: is this the type that we want?
                fwdpy11::discrete_demography::DemographyError);
        }
    // pop hasn't evolved!
    BOOST_REQUIRE_EQUAL(pop.generation, 0);
}

BOOST_FIXTURE_TEST_CASE(test_initial_pop_size_invalid_island_model,
                        common_setup<TwoDemePerpetualIslandModel>)
{
    {
        // The demes model specifies N = [100, 100].
        // Here, we will start with N for one deme
        // This is a hard error!
        // We test values ~100 because TDD found cases of memory errors
        // when N >> the correct N but things can pass when N =~ the
        // correct N. Gotta love UB :).
        std::vector<std::uint32_t> initial_popsizes{50, 99, 101, 200};
        for (auto initial_n : initial_popsizes)
            {
                // reset the fixture
                pop = fwdpy11::DiploidPopulation(initial_n, 10.0);
                auto model = build_model();
                fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

                BOOST_CHECK_THROW(
                    {
                        evolve_with_tree_sequences(
                            rng, pop, recorder, 10, forward_demes_graph, 10, 0., 0.,
                            mregions, recregions, gvalue_ptrs, sample_recorder_callback,
                            stopping_criterion, post_simplification_recorder, options);
                    },
                    // TODO: is this the type that we want?
                    fwdpy11::discrete_demography::DemographyError);
            }
        // pop hasn't evolved!
        BOOST_REQUIRE_EQUAL(pop.generation, 0);
    }

    {
        std::vector<std::vector<std::uint32_t>> initial_popsizes{
            {100, 0}, {0, 100}, {10, 100}, {100, 99, 1}};
        for (auto initial_n : initial_popsizes)
            {
                // reset the fixture
                pop = fwdpy11::DiploidPopulation(initial_n, 10.0);
                auto model = build_model();
                fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
                BOOST_CHECK_THROW(
                    {
                        evolve_with_tree_sequences(
                            rng, pop, recorder, 10, forward_demes_graph, 10, 0., 0.,
                            mregions, recregions, gvalue_ptrs, sample_recorder_callback,
                            stopping_criterion, post_simplification_recorder, options);
                    },
                    // TODO: is this the type that we want?
                    fwdpy11::discrete_demography::DemographyError);
            }
        // pop hasn't evolved!
        BOOST_REQUIRE_EQUAL(pop.generation, 0);
    }
}

BOOST_FIXTURE_TEST_CASE(test_generation_time_past_end_of_model,
                        common_setup<SingleDemeModel>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // The model is only set to simulate 10 generations
    pop.generation = 11;
    BOOST_CHECK_THROW(
        {
            evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 10,
                                       0., 0., mregions, recregions, gvalue_ptrs,
                                       sample_recorder_callback, stopping_criterion,
                                       post_simplification_recorder, options);
        },
        fwdpy11::discrete_demography::DemographyError);
    BOOST_REQUIRE_EQUAL(pop.generation, 11);
}

BOOST_FIXTURE_TEST_CASE(
    test_size_history_single_deme_model_one_size_change,
    common_setup_with_sample_history_recording<SingleDemeModelOneSizeChange>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 60);
    BOOST_REQUIRE_EQUAL(size_history.size(), 1);
    for (auto& x : size_history)
        {
            BOOST_REQUIRE_EQUAL(x.second.size(), 60);
            auto ngens_at_100 = std::count_if(std::begin(x.second), std::end(x.second),
                                              [](auto& x) { return x.second == 100; });
            BOOST_REQUIRE_EQUAL(ngens_at_100, 10);
            BOOST_REQUIRE_EQUAL(std::count_if(std::begin(x.second), std::end(x.second),
                                              [](auto& x) { return x.second == 200; }),
                                50);
            BOOST_REQUIRE(std::all_of(std::begin(x.second),
                                      std::begin(x.second) + ngens_at_100,
                                      [](auto& y) { return y.second == 100; }));
            BOOST_REQUIRE(std::all_of(std::begin(x.second) + ngens_at_100,
                                      std::end(x.second),
                                      [](auto& y) { return y.second == 200; }));
        }
}

BOOST_FIXTURE_TEST_CASE(
    test_size_history_two_deme_island_model_with_size_change_and_excinction,
    common_setup_with_sample_history_recording<
        TwoDemePerpetualIslandModelWithSizeChangeAndExtinction>)
{
    pop = fwdpy11::DiploidPopulation({100, 100}, 10.);
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 60);
    BOOST_REQUIRE_EQUAL(size_history.size(), 2);
    BOOST_REQUIRE(std::all_of(std::begin(pop.diploid_metadata),
                              std::end(pop.diploid_metadata),
                              [](const auto& md) { return md.deme == 1; }));

    // Deme 0 went exctinct 25 gens ago
    BOOST_REQUIRE_EQUAL(size_history[0].size(), 60 - 25);
    BOOST_REQUIRE_EQUAL(size_history[1].size(), 60);
    BOOST_REQUIRE(std::all_of(std::begin(size_history[0]), std::end(size_history[0]),
                              [](const auto& y) { return y.second == 100; }));

    // Deme 1 changed from N=100 to N=200 50 gens ago
    BOOST_REQUIRE(std::all_of(std::begin(size_history[1]),
                              std::begin(size_history[1]) + 60 - 50,
                              [](const auto& y) { return y.second == 100; }));
    BOOST_REQUIRE(std::all_of(std::begin(size_history[1]) + 60 - 50,
                              std::end(size_history[1]),
                              [](const auto& y) { return y.second == 200; }));
}

BOOST_FIXTURE_TEST_CASE(test_size_history_two_demes_unequal_merge,
                        common_setup_with_sample_history_recording<TwoDemesUnequalMerge>)
{
    pop = fwdpy11::DiploidPopulation({100, 75}, 10.);
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 35);
    BOOST_REQUIRE_EQUAL(size_history.size(), 3);

    BOOST_REQUIRE(std::all_of(std::begin(pop.diploid_metadata),
                              std::end(pop.diploid_metadata),
                              [](const auto& md) { return md.deme == 2; }));

    // Deme 0 went extinct 25 gens ago
    BOOST_REQUIRE_EQUAL(size_history[0].size(), 35 - 25);
    BOOST_REQUIRE(std::all_of(std::begin(size_history[0]), std::end(size_history[0]),
                              [](const auto& i) { return i.first < 35 - 25 + 1; }));
    // Deme 1 went extinct 10 gens ago
    BOOST_REQUIRE_EQUAL(size_history[1].size(), 35 - 10);
    BOOST_REQUIRE(std::all_of(std::begin(size_history[1]), std::end(size_history[1]),
                              [](const auto& i) { return i.first < 35 - 10 + 1; }));
    // Deme 2 first appeared 25 gens ago
    BOOST_REQUIRE_EQUAL(size_history[2].size(), 25);
    BOOST_REQUIRE(std::all_of(std::begin(size_history[2]), std::end(size_history[2]),
                              [](const auto& i) { return i.first >= 35 - 25 + 1; }));
}

// NOTE: upstream has more extensive tests of correctness
BOOST_FIXTURE_TEST_CASE(test_linear_size_change_history,
                        common_setup_with_sample_history_recording<LinearSizeChange>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 30);
    BOOST_REQUIRE_EQUAL(pop.N, 55);

    BOOST_REQUIRE_EQUAL(size_history[0][9].second, 100);
    // this is wrong: need to calc slope, etc..
    double slope = (200. - 100.) / 10.;
    BOOST_REQUIRE_EQUAL(size_history[0][10].second, static_cast<unsigned>(100. + slope));
    BOOST_REQUIRE_EQUAL(size_history[0][15].second,
                        static_cast<unsigned>(100. + 6. * slope));
    BOOST_REQUIRE_EQUAL(size_history[0][18].first, 19);
    BOOST_REQUIRE_EQUAL(size_history[0][18].second,
                        static_cast<unsigned>(100. + 9. * slope));
    BOOST_REQUIRE_EQUAL(size_history[0][19].first, 20);
    BOOST_REQUIRE_EQUAL(size_history[0][19].second,
                        static_cast<unsigned>(100. + 10. * slope));
    // Pop size is flat at post-growth size
    // for rest of model
    for (unsigned i = 21; i < 31; ++i)
        {
            BOOST_REQUIRE_EQUAL(size_history[0][i - 1].first, i);
            BOOST_REQUIRE_EQUAL(size_history[0][i - 1].second, 55);
        }
}

BOOST_FIXTURE_TEST_CASE(
    test_ancestry_proportions_from_pulse_with_burnin,
    common_setup_with_ancestry_tracking<VeryRecentPulseTwoGenerationsAgo>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    BOOST_REQUIRE_EQUAL(forward_demes_graph.number_of_demes(), 2);
    pop = fwdpy11::DiploidPopulation({50, 50}, 1.);
    update_parents(pop);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 12);

    auto condition = [](std::uint32_t a) { return a == 50; };

    // pre-pulse
    for (std::uint32_t g = 1; g < 10; ++g)
        {
            auto found = validate_ancestry(this->ancestry, g, 0, 0, condition);
            BOOST_REQUIRE(found);
            found = validate_ancestry(this->ancestry, g, 1, 1, condition);
            BOOST_REQUIRE(found);
        }
    auto found = validate_ancestry(this->ancestry, 10, 0, 1, condition);
    BOOST_REQUIRE(found);
    found = validate_ancestry(this->ancestry, 10, 1, 0, condition);
    BOOST_REQUIRE(found);
    found = validate_ancestry(this->ancestry, 11, 0, 0, condition);
    BOOST_REQUIRE(found);
    found = validate_ancestry(this->ancestry, 11, 1, 1, condition);
    BOOST_REQUIRE(found);
}

BOOST_FIXTURE_TEST_CASE(
    test_ancestry_with_extreme_migration_until_one_generation_ago,
    common_setup_with_ancestry_tracking<ExtremeMigrationUntilOneGenerationAgo>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    BOOST_REQUIRE_EQUAL(forward_demes_graph.number_of_demes(), 2);
    pop = fwdpy11::DiploidPopulation({1000, 1000}, 1.);
    update_parents(pop);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 60, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);

    auto validate_complete_ancestry = [](std::uint32_t p) { return p == 1000; };

    auto validate_partial_ancestry = [](std::uint32_t p) { return p > 0 && p < 1000; };

    // NOTE: at the deme sizes in this model
    // and migration rates reflecting 50% migrant
    // ancestry each generation, is is very unlikely
    // to not observe some ancestry each generation
    // in each combo tested below.
    for (std::uint32_t g = 1; g < 10; ++g)
        {
            for (std::int32_t i = 0; i < 2; ++i)
                {
                    for (std::int32_t j = 0; j < 2; ++j)
                        {
                            auto found = validate_ancestry(ancestry, g, i, j,
                                                           validate_partial_ancestry);
                            BOOST_REQUIRE(found);
                        }
                }
        }
    auto found = validate_ancestry(ancestry, 10, 0, 0, validate_complete_ancestry);
    BOOST_REQUIRE(found);

    found = validate_ancestry(ancestry, 10, 1, 1, validate_complete_ancestry);
    BOOST_REQUIRE(found);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_evolvets_start_stop_start)

// Tests of repeatedly entering in to the same model
BOOST_FIXTURE_TEST_CASE(
    test_ancestry_with_extreme_migration_until_one_generation_ago_start_stop,
    common_setup_with_ancestry_tracking<ExtremeMigrationUntilOneGenerationAgo>)
{
    // NOTE: this is a copy-paste of a test from above.
    // The only difference is that we call into the
    // evolve function for one generation at a time.
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    BOOST_REQUIRE_EQUAL(forward_demes_graph.number_of_demes(), 2);
    pop = fwdpy11::DiploidPopulation({1000, 1000}, 1.);
    update_parents(pop);
    auto last_generation = pop.generation;
    for (std::size_t i = 0; i < 10; ++i)
        {
            evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 1,
                                       0., 0., mregions, recregions, gvalue_ptrs,
                                       sample_recorder_callback, stopping_criterion,
                                       post_simplification_recorder, options);
            BOOST_REQUIRE_EQUAL(pop.generation, last_generation + 1);
            last_generation += 1;
        }
    BOOST_REQUIRE_EQUAL(pop.generation, 10);

    auto validate_complete_ancestry = [](std::uint32_t p) { return p == 1000; };

    auto validate_partial_ancestry = [](std::uint32_t p) { return p > 0 && p < 1000; };

    // NOTE: at the deme sizes in this model
    // and migration rates reflecting 50% migrant
    // ancestry each generation, is is very unlikely
    // to not observe some ancestry each generation
    // in each combo tested below.
    for (std::uint32_t g = 1; g < 10; ++g)
        {
            for (std::int32_t i = 0; i < 2; ++i)
                {
                    for (std::int32_t j = 0; j < 2; ++j)
                        {
                            auto found = validate_ancestry(ancestry, g, i, j,
                                                           validate_partial_ancestry);
                            BOOST_REQUIRE(found);
                        }
                }
        }
    auto found = validate_ancestry(ancestry, 10, 0, 0, validate_complete_ancestry);
    BOOST_REQUIRE(found);

    found = validate_ancestry(ancestry, 10, 1, 1, validate_complete_ancestry);
    BOOST_REQUIRE(found);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_residual_selfing)

BOOST_FIXTURE_TEST_CASE(test_no_residual_selfing, common_setup<DemeSizeIsOne>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    pop = fwdpy11::DiploidPopulation(1, 1.0);
    evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 1, 0., 0.,
                               mregions, recregions, gvalue_ptrs,
                               sample_recorder_callback, stopping_criterion,
                               post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 1);
}

BOOST_FIXTURE_TEST_CASE(test_no_residual_selfing_exception, common_setup<DemeSizeIsOne>)
{
    auto model = build_model();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);
    pop = fwdpy11::DiploidPopulation(1, 1.0);
    options.allow_residual_selfing = false;
    BOOST_REQUIRE_THROW(
        {
            evolve_with_tree_sequences(rng, pop, recorder, 10, forward_demes_graph, 1,
                                       0., 0., mregions, recregions, gvalue_ptrs,
                                       sample_recorder_callback, stopping_criterion,
                                       post_simplification_recorder, options);
        },
        fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_SUITE_END()
