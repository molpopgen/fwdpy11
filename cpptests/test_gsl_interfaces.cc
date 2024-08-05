#include <boost/test/unit_test.hpp>

#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
#include <core/gsl/gsl_discrete.hpp>
#include <fwdpy11/rng.hpp>
#include <gsl/gsl_randist.h>
#include <stdexcept>

BOOST_AUTO_TEST_SUITE(test_gsl_interfaces)

BOOST_AUTO_TEST_CASE(test_negative_values_in_lookup_table)
{
    fwdpy11::gsl_scoped_convert_error_to_exception h;
    std::vector<double> weights{-1., 1., 1.};
    BOOST_REQUIRE_EQUAL(weights.size(), 3);
    fwdpp::gsl_ran_discrete_t_ptr lookup;

    BOOST_REQUIRE_THROW(
        { fwdpy11_core::update_lookup_table(weights.data(), weights.size(), lookup); },
        fwdpy11::GSLError);
}

BOOST_AUTO_TEST_CASE(test_negative_values_in_lookup_table_manual_catch)
{
    fwdpy11::gsl_scoped_disable_error_handler_wrapper h;
    std::vector<double> weights{-1., 1., 1.};
    BOOST_REQUIRE_EQUAL(weights.size(), 3);
    fwdpp::gsl_ran_discrete_t_ptr lookup;

    BOOST_REQUIRE_THROW(
        { fwdpy11_core::update_lookup_table(weights.data(), weights.size(), lookup); },
        fwdpy11::GSLError);
}

// This tests a dangerous edge case: if all fitnesses are 0,
// the simulation is effectively neutral b/c all samples have the same weight!
BOOST_AUTO_TEST_CASE(test_all_zeros_in_lookup_table)
{
    fwdpy11::gsl_scoped_disable_error_handler_wrapper h;
    std::vector<double> weights{0., 0., 0.};
    BOOST_REQUIRE_EQUAL(weights.size(), 3);

    auto table = gsl_ran_discrete_preproc(weights.size(), weights.data());
    if (table == nullptr)
        {
            throw std::runtime_error("null table due to all zeros");
        }
    fwdpp::gsl_ran_discrete_t_ptr lookup(table);
}

BOOST_AUTO_TEST_CASE(test_catch_all_zeros_in_lookup_table)
{
    fwdpy11::gsl_scoped_disable_error_handler_wrapper h;
    std::vector<double> weights{0., 0., 0.};
    BOOST_REQUIRE_EQUAL(weights.size(), 3);
    fwdpp::gsl_ran_discrete_t_ptr lookup;

    BOOST_REQUIRE_THROW(
        { fwdpy11_core::update_lookup_table(weights.data(), weights.size(), lookup); },
        fwdpy11::GSLError);
}

BOOST_AUTO_TEST_SUITE_END()
