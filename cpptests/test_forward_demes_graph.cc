#include <boost/test/unit_test.hpp>
#include <fwdpy11/discrete_demography/forward_demes_graph.hpp>

BOOST_AUTO_TEST_SUITE(test_forward_demes_graph)

BOOST_AUTO_TEST_CASE(create_empty_graph)
{
    fwdpy11::discrete_demography::ForwardDemesGraph g;
}

BOOST_AUTO_TEST_SUITE_END()
