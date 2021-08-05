#include "fwdpy11/discrete_demography/MigrationMatrix.hpp"
#include "fwdpy11/discrete_demography/SetDemeSize.hpp"
#include "fwdpy11/discrete_demography/SetExponentialGrowth.hpp"
#include "fwdpy11/discrete_demography/SetMigrationRates.hpp"
#include "fwdpy11/discrete_demography/SetSelfingRate.hpp"
#include <boost/test/unit_test.hpp>
#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <stdexcept>
#include <vector>

#include "discrete_demography_util.hpp"
#include "fwdpy11/discrete_demography/constants.hpp"
#include "discrete_demography_fixtures.hpp"

struct MigrationMatrix_fixture
{
    fwdpy11::discrete_demography::MigrationMatrix migmatrix;

    MigrationMatrix_fixture() : migmatrix{std::vector<double>{1, 0, 0, 1}, 2, false}
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_lowlevel_discrete_demography_objects)

BOOST_AUTO_TEST_CASE(test_SetDemeSize)
{
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetDemeSize(0, -1, 10, false); },
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_SetExponentialGrowth)
{
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetExponentialGrowth(0, -1, 1.); },
                      std::invalid_argument);
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetExponentialGrowth(0, 1, -1.); },
                      std::invalid_argument);
    BOOST_CHECK_THROW(
        { fwdpy11::discrete_demography::SetExponentialGrowth(0, 1, 1. / 0.); },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_SetSelfingRate)
{
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetSelfingRate(0, -1, 1.); },
                      std::invalid_argument);
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetSelfingRate(0, 1, 1.1); },
                      std::invalid_argument);
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetSelfingRate(0, 1, -1.); },
                      std::invalid_argument);
    BOOST_CHECK_THROW({ fwdpy11::discrete_demography::SetSelfingRate(0, 1, 1. / 0.); },
                      std::invalid_argument);
}

// NOTE: this only tests the minimal interface
BOOST_AUTO_TEST_CASE(test_MassMigrations)
{
    BOOST_CHECK_THROW({ move_individuals(0, 0, 1, 0. / 0., true); },
                      std::invalid_argument);
    BOOST_CHECK_THROW({ move_individuals(0, 0, 1, 2., true); }, std::invalid_argument);
    BOOST_CHECK_THROW({ move_individuals(0, 0, 1, -1., true); }, std::invalid_argument);
    BOOST_CHECK_THROW({ move_individuals(0, -1, 1, 1., true); }, std::invalid_argument);
    BOOST_CHECK_THROW({ move_individuals(0, 0, -1, 1., true); }, std::invalid_argument);
    BOOST_CHECK_THROW({ move_individuals(0, 0, 0, 1., true); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_MigrationMatrix)
{
    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::MigrationMatrix(std::vector<double>{0, 1, 0},
                                                          2, false);
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::MigrationMatrix(
                std::vector<double>{0, 1, 0, -1}, 2, false);
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::MigrationMatrix(
                std::vector<double>{0, 1, 1, 1}, 2, false);
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::MigrationMatrix(
                std::vector<double>{0, 1, 0, 1. / 0.}, 2, false);
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_SetMigrationRates)
{
    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(0, -1,
                                                            std::vector<double>{1, 0});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(0, 1,
                                                            std::vector<double>{1, 1});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(0, 1, std::vector<double>{});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, 1, std::vector<double>{0.1, 0.25});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, 1, std::vector<double>{0.1, 0.25 / 0.});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(0, 1,
                                                            std::vector<double>{-1, 1});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, std::vector<double>{0.1, 0.25});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        { fwdpy11::discrete_demography::SetMigrationRates(0, std::vector<double>{}); },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, std::vector<double>{0.1, 0.9, 0.1, 0. / 0.});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, std::vector<double>{0.1, 0.9, 0.1, -0.1});
        },
        std::invalid_argument);

    BOOST_CHECK_THROW(
        {
            fwdpy11::discrete_demography::SetMigrationRates(
                0, std::vector<double>{0.1, 0.9, 0.1, 0.1});
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_set_migration_rates, MigrationMatrix_fixture)

BOOST_AUTO_TEST_CASE(set_row)
{
    BOOST_CHECK_NO_THROW(try {
        migmatrix.set_migration_rates(0, std::vector<double>{0.5, 0.5});
    } catch (...){});
}

BOOST_AUTO_TEST_CASE(set_matrix)
{
    BOOST_CHECK_NO_THROW(try {
        migmatrix.set_migration_rates(fwdpy11::discrete_demography::NULLDEME,
                                      std::vector<double>{0.5, 0.5, 0.5, 0.5});
    } catch (...){});
}

BOOST_AUTO_TEST_CASE(bad_set_matrix)
{
    BOOST_CHECK_THROW(
        {
            migmatrix.set_migration_rates(fwdpy11::discrete_demography::NULLDEME,
                                          std::vector<double>{0.5, 0.5});
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(source_deme_id_out_of_range)
{
    BOOST_CHECK_THROW(
        {
            migmatrix.set_migration_rates(
                10, std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1});
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(wrong_number_of_rates)
{
    BOOST_CHECK_THROW(
        {
            migmatrix.set_migration_rates(
                0, std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1});
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_DiscreteDemography, population_fixture)

BOOST_AUTO_TEST_CASE(set_deme_size_twice)
{
    set_deme_sizes.emplace_back(0, 1, 100, true);
    set_deme_sizes.emplace_back(0, 1, 10, false);

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_growth_rate_twice)
{
    set_growth_rates.emplace_back(0, 1, 1.1);
    set_growth_rates.emplace_back(0, 1, 1.2);

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_selfing_rates_twice)
{
    set_selfing_rates.emplace_back(0, 1, 0.1);
    set_selfing_rates.emplace_back(0, 1, 0.2);

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_migration_rates_twice)
{
    set_migration_rates.emplace_back(0, 1, std::vector<double>{0.5, 0.5});
    set_migration_rates.emplace_back(0, 1, std::vector<double>{0.25, 0.75});

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(set_migration_rates_but_there_is_no_migration_matrix)
{
    set_migration_rates.emplace_back(0, 1, std::vector<double>{0.5, 0.5});

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(wrong_size_migration_matrix)
{
    auto mptr = make_gsl_matrix(3);
    gsl_matrix_set_identity(mptr.get());
    set_migmatrix(convert_matrix(mptr), 3, false);
    set_migration_rates.emplace_back(fwdpy11::discrete_demography::NULLDEME,
                                     std::vector<double>{0.5, 0.5, 0.5, 0.5});

    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(move_too_many_individuals)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.25, true));
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.85, true));
    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(mass_migrate_twice)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.25, true));
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.35, true));
    BOOST_CHECK_THROW({ auto ddemog = make_model(); }, std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

