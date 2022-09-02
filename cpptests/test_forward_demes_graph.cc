#include <boost/test/unit_test.hpp>

#include <boost/test/unit_test_suite.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <core/demes/forward_graph.hpp>
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
    g.initialize_model(pop.generation);
    auto end_time = g.model_end_time();
    BOOST_REQUIRE_EQUAL(end_time, 11);
    while (g.iterating_model())
        {
            g.iterate_state();
            pop.generation += 1;
        }
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_FIXTURE_TEST_CASE(single_deme_model_with_burn_in_THE_CODE_WE_WANT_TO_WRITE,
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
            g.iterate_state();
            pop.generation += 1;
        }
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_AUTO_TEST_SUITE_END()
