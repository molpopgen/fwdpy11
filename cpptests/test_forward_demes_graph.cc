#include <boost/test/unit_test.hpp>

#include <boost/test/unit_test_suite.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <core/demes/forward_graph.hpp>
#include <numeric>
#include <stdexcept>

#include "forward_demes_graph_fixtures.hpp"

BOOST_AUTO_TEST_SUITE(test_forward_demes_graph)

BOOST_FIXTURE_TEST_CASE(test_zero_length_single_deme_model, SingleDemeModel)
{
    fwdpy11_core::ForwardDemesGraph g(yaml, 0);
    fwdpy11::DiploidPopulation pop(100, 1.0);
    g.initialize_model(pop.generation);
    auto end_time = g.model_end_time();
    BOOST_REQUIRE_EQUAL(end_time, 1);
    while (g.iterating_model())
        {
            g.iterate_state();
            pop.generation += 1;
        }
    BOOST_REQUIRE_EQUAL(pop.generation, 0);
}

BOOST_FIXTURE_TEST_CASE(single_deme_model_with_burn_in, SingleDemeModel)
{
    fwdpy11_core::ForwardDemesGraph g(yaml, 10);
    BOOST_REQUIRE_EQUAL(g.number_of_demes(), 1);
    fwdpy11::DiploidPopulation pop(100, 1.0);
    BOOST_CHECK_EQUAL(pop.N, g.sum_deme_sizes_at_time_zero());
    g.initialize_model(pop.generation);
    auto end_time = g.model_end_time();
    BOOST_REQUIRE_EQUAL(end_time, 11);
    while (g.iterating_model())
        {
            auto num = 0;
            for (auto p : g.parental_deme_sizes())
                {
                    BOOST_REQUIRE_EQUAL(p, 100.0);
                    num += 1;
                }
            BOOST_REQUIRE_EQUAL(num, g.number_of_demes());
            num = 0;
            for (auto p : g.offspring_deme_sizes())
                {
                    BOOST_REQUIRE_EQUAL(p, 100.0);
                    if (p > 0.0)
                        {
                            auto ancestry = g.offspring_ancestry_proportions(num);
                            auto sum_proportions = std::accumulate(
                                std::begin(ancestry), std::end(ancestry), 0.0);
                            BOOST_REQUIRE(std::fabs(sum_proportions - 1.0)
                                          <= std::numeric_limits<double>::epsilon());
                        }
                    num += 1;
                }
            BOOST_REQUIRE_EQUAL(num, g.number_of_demes());
            num = 0;
            for (auto s : g.offspring_selfing_rates())
                {
                    BOOST_REQUIRE_EQUAL(s, 0.0);
                    num += 1;
                }
            BOOST_REQUIRE_EQUAL(num, g.number_of_demes());
            num = 0;
            for (auto c : g.offspring_cloning_rates())
                {
                    BOOST_REQUIRE_EQUAL(c, 0.0);
                    num += 1;
                }
            BOOST_REQUIRE_EQUAL(num, g.number_of_demes());

            pop.generation += 1;
            g.iterate_state();
        }
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_FIXTURE_TEST_CASE(single_deme_model_with_burn_in_invalid_ancestry_request,
                        SingleDemeModel)
{
    fwdpy11_core::ForwardDemesGraph g(yaml, 10);
    BOOST_REQUIRE_EQUAL(g.number_of_demes(), 1);
    fwdpy11::DiploidPopulation pop(100, 1.0);
    g.initialize_model(pop.generation);
    auto end_time = g.model_end_time();
    BOOST_REQUIRE_EQUAL(end_time, 11);
    while (g.iterating_model())
        {
            // this is an invalid index
            BOOST_REQUIRE_THROW(
                {
                    // We are expecting (?) the API to set the error code here.
                    /* auto ancestry = */ g.offspring_ancestry_proportions(
                        g.number_of_demes());
                },
                fwdpy11::discrete_demography::DemographyError);
            g.iterate_state();
        }
}

BOOST_FIXTURE_TEST_CASE(has_non_integer_start_size, NonIntegerStartSize)
{
    BOOST_REQUIRE_THROW({ fwdpy11_core::ForwardDemesGraph g(yaml, 10); },
                        std::invalid_argument);
}

BOOST_FIXTURE_TEST_CASE(has_non_integer_end_size, NonIntegerEndSize)
{
    BOOST_REQUIRE_THROW({ fwdpy11_core::ForwardDemesGraph g(yaml, 10); },
                        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_forward_demes_graph_with_bad_models)

BOOST_FIXTURE_TEST_CASE(test_bad_epoch_rounding_02, BadEpochRounding02)
{
    BOOST_REQUIRE_THROW({ fwdpy11_core::ForwardDemesGraph g(yaml, 10); },
                        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
