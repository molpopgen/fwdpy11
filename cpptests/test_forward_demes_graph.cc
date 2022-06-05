#include <stdexcept>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <fwdpy11/discrete_demography/forward_demes_graph.hpp>

BOOST_AUTO_TEST_SUITE(test_forward_demes_graph)

BOOST_AUTO_TEST_CASE(add_deme)
{
    fwdpy11::discrete_demography::ForwardDemesGraph g;

    // return pointer/reference to deme with a given start time.
    // arguments are name, id, start_time
    auto deme = g.add_deme("CEU", 0, 0);
}

BOOST_AUTO_TEST_CASE(add_deme_that_already_exists)
{
    fwdpy11::discrete_demography::ForwardDemesGraph g;

    g.add_deme("CEU", 0, 0);
    BOOST_CHECK_THROW({ g.add_deme("CEU", 0, 0); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(add_epoch_to_deme)
{
    fwdpy11::discrete_demography::ForwardDemesGraph g;

    auto deme = g.add_deme("CEU", 0, 0);

    // end tie, start size, end size, cloning, selfing, size_function
    deme.add_epoch(100, 150, 100, 0.,
                   fwdpy11::discrete_demography::WrightFisherSelfing(),
                   fwdpy11::discrete_demography::ConstantSizeFunction());
}

BOOST_AUTO_TEST_SUITE_END()
