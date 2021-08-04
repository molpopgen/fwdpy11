#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <vector>
#include "discrete_demography_fixtures.hpp"
#include "discrete_demography_util.hpp"
#include "fwdpy11/discrete_demography/exceptions.hpp"

BOOST_FIXTURE_TEST_SUITE(test_lowlevel_discrete_demography_bad_models,
                         population_fixture)

BOOST_AUTO_TEST_CASE(test_migration_matrix_too_small)
/*
 * The mig matrix is 2x2 but a mass migration event
 * tries to create a 3rd deme, which triggers an exception
 */
{
    set_migmatrix(std::vector<double>{0.5, 0.5, 0.5, 0.5}, 2, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 5); },
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_migration_matrix_too_large)
// One deme and a 2x2 matrix
{
    set_migmatrix(std::vector<double>{0, 1, 1, 0}, 2, false);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 5); },
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_change_migration_rates_simple_two_deme_migration_bad_matrix)
/*
 * For a 2-deme model, the mig matrix is
 * [0, 1
 *  1, 0]
 * so that all offspring have both parents from the other deme,
 * which gives us an easy check on how many migration events
 * will be recorded by the test simulation.
 * 
 * After 3 generations, we reset the migration rates to be
 * [[0.5, 0.5],
 *  [0, 0]],
 * which leads to there being no parents for deme 1, raising a
 * fwdpy11.DemographyError exception.
 */
{
    set_migmatrix(std::vector<double>{0, 1, 1, 0}, 2, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    set_migration_rates.emplace_back(3, std::vector<double>{0.5, 0.5, 0, 0});
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 5); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_growth_in_deme_that_doesnt_exist)
{
    set_growth_rates.emplace_back(0, 1, 0.5);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_growth_in_deme_that_went_extinct)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    set_deme_sizes.emplace_back(5, 1, 0, true);
    set_growth_rates.emplace_back(6, 1, 0.5);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_set_selfing_in_deme_that_doesnt_exist)
{
    set_selfing_rates.emplace_back(0, 1, 0.5);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_set_selfing_in_deme_that_went_extinct)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    set_deme_sizes.emplace_back(5, 1, 0, true);
    set_selfing_rates.emplace_back(6, 1, 0.5);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_migration_into_empty_deme)
{
    auto mptr = make_gsl_matrix(3);
    auto cview = gsl_matrix_column(mptr.get(), 2);
    for (std::size_t i = 0; i < cview.vector.size; ++i)
        {
            gsl_vector_set(&cview.vector, i, 1.0);
        }
    auto v = convert_matrix(mptr);
    set_migmatrix(convert_matrix(mptr), 3, false);

    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));

    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_migration_from_an_empty_deme)
{
    auto mptr = make_gsl_matrix(3);
    std::fill(mptr->data, mptr->data + mptr->size1 * mptr->size2, 0.25);
    for (std::size_t i = 0; i < 3; ++i)
        {
            gsl_matrix_set(mptr.get(), i, i, 0.5);
        }
    set_migmatrix(convert_matrix(mptr), 3, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_migration_from_an_empty_deme_into_an_empty_deme)
{
    auto mptr = make_gsl_matrix(3);
    gsl_matrix_set(mptr.get(), 0, 1, 1.0);
    gsl_matrix_set(mptr.get(), 2, 2, 1.0);
    set_migmatrix(convert_matrix(mptr), 3, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_no_valid_parents)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 1, true));
    set_deme_sizes.emplace_back(0, 0, 100, true);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_no_valid_parents_alternate_method)
{
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 1, true));
    set_deme_sizes.emplace_back(0, 0, 0, true);
    set_deme_sizes.emplace_back(2, 0, 100, true);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_no_valid_parents_with_migration)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 1, true));
    set_deme_sizes.emplace_back(0, 0, 100, true);
    set_migmatrix(std::vector<double>{0.5, 0.5, 0.5, 0.5}, 2, false);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(test_no_valid_parents_with_migration_v2)
/*
 * This model works because we take pains to 
 * ensure that there are possible parents for all demes
 * at all times
 */
{
    // Set up initial migration: equal rates all over!
    set_migmatrix(std::vector<double>{0.5, 0.5, 0.5, 0.5}, 2, false);
    // Move everyone immediately -- the migmatrix is now invalid!
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 1, true));
    // But we step in and set deme 0's size back to 100
    set_deme_sizes.emplace_back(0, 0, 100, true);
    // And, we make sure that deme 1 provides 100% of the ancestry to deme 0
    // and to deme 1.
    set_migration_rates.emplace_back(0, std::vector<double>{0, 1, 0, 1});
    auto ddemog = make_model();
    DiscreteDemography_roundtrip(rng, pop, ddemog, 20);
    BOOST_CHECK_NO_THROW(
        try { DiscreteDemography_roundtrip(rng, pop, ddemog, 20); } catch (...){});
}

BOOST_AUTO_TEST_CASE(move_from_empty_deme)
{
    set_deme_sizes.emplace_back(0, 0, 0, true);
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 1., true));
    mass_migrations.emplace_back(move_individuals(1, 0, 1, 0.1, true));
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_CASE(copy_from_empty_deme)
{
    set_deme_sizes.emplace_back(0, 0, 0, true);
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 1., true));
    mass_migrations.emplace_back(copy_individuals(1, 0, 1, 0.1, true));
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::DemographyError);
}

BOOST_AUTO_TEST_SUITE_END()
