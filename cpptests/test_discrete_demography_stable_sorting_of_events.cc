#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <boost/test/unit_test.hpp>

struct demographic_event_stable_sorting_fixture
{
    fwdpy11::discrete_demography::DiscreteDemography demog;

    fwdpy11::discrete_demography::DiscreteDemography
    init()
    // This function does not create a valid model that can be simulated.
    {
        std::vector<fwdpy11::discrete_demography::MassMigration> massmigs;
        std::vector<fwdpy11::discrete_demography::SetExponentialGrowth> growth;
        std::vector<fwdpy11::discrete_demography::SetDemeSize> sizechanges;
        std::vector<fwdpy11::discrete_demography::SetSelfingRate> selfing;
        std::vector<fwdpy11::discrete_demography::SetMigrationRates> migration;

        massmigs.emplace_back(5, 5, 3, 0, -1, 0.25, true, false, true);
        massmigs.emplace_back(5, 3, 2, 0, -1, 0.25, true, false, true);
        massmigs.emplace_back(4, 4, 3, 0, -1, 0.25, true, false, true);
        massmigs.emplace_back(4, 1, 2, 0, -1, 0.25, true, false, true);

        growth.emplace_back(10, 1, 0.5);
        growth.emplace_back(10, 0, 0.25);
        growth.emplace_back(5, 0, 1.);
        growth.emplace_back(5, 6, 0.1);
        growth.emplace_back(5, 3, 0.2);

        sizechanges.emplace_back(10, 1, 11, true);
        sizechanges.emplace_back(10, 0, 300, true);
        sizechanges.emplace_back(5, 0, 33, true);
        sizechanges.emplace_back(5, 6, 1, true);
        sizechanges.emplace_back(5, 3, 22, true);

        selfing.emplace_back(10, 1, 0.5);
        selfing.emplace_back(10, 0, 0.25);
        selfing.emplace_back(5, 0, 1.);
        selfing.emplace_back(5, 6, 0.);
        selfing.emplace_back(5, 3, 0.2);

        std::vector<double> migmatrix(5 * 5, 0.);
        for (std::size_t i = 0; i < 5; ++i)
            {
                migmatrix[5 * i + i] = 1.;
            }

        migration.emplace_back(5, 3, std::vector<double>{0., 0., 1., 0., 0.});
        migration.emplace_back(4, 4, std::vector<double>{0., 0., 1., 0., 0.});
        migration.emplace_back(4, migmatrix);

        return fwdpy11::discrete_demography::DiscreteDemography(
            std::move(massmigs), std::move(growth), std::move(sizechanges),
            std::move(selfing),
            fwdpy11::discrete_demography::MigrationMatrix(std::move(migmatrix), 5,
                                                          false),
            std::move(migration));
    }

    demographic_event_stable_sorting_fixture() : demog{init()}
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(test_stable_sorting_of_demographic_events,
                         demographic_event_stable_sorting_fixture)

BOOST_AUTO_TEST_CASE(test_mass_migrations)
{
    std::vector<std::uint32_t> when{4, 4, 5, 5};
    std::vector<int> source{4, 1, 5, 3};
    std::vector<int> dest{3, 2, 3, 2};
    for (std::size_t i = 0; i < demog.get_mass_migrations().size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(demog.get_mass_migrations()[i].when, when[i]);
            BOOST_REQUIRE_EQUAL(demog.get_mass_migrations()[i].source, source[i]);
            BOOST_REQUIRE_EQUAL(demog.get_mass_migrations()[i].destination, dest[i]);
        }
}

BOOST_AUTO_TEST_CASE(test_set_migration_rates)
{
    std::vector<std::uint32_t> when{4, 4, 5};
    std::vector<int> deme{4, fwdpy11::discrete_demography::NULLDEME, 3};
    for (std::size_t i = 0; i < demog.get_set_migration_rates().size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(demog.get_set_migration_rates()[i].when, when[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_migration_rates()[i].deme, deme[i]);
        }
}

BOOST_AUTO_TEST_CASE(test_size_changes)
{
    std::vector<std::uint32_t> when{5, 5, 5, 10, 10};
    std::vector<int> deme{0, 6, 3, 1, 0};
    std::vector<std::uint32_t> new_size{33, 1, 22, 11, 300};

    for (std::size_t i = 0; i < demog.get_set_deme_sizes().size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(demog.get_set_deme_sizes()[i].when, when[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_deme_sizes()[i].deme, deme[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_deme_sizes()[i].new_size, new_size[i]);
        }
}

BOOST_AUTO_TEST_CASE(test_growth_rate_changes)
{
    std::vector<std::uint32_t> when{5, 5, 5, 10, 10};
    std::vector<int> deme{0, 6, 3, 1, 0};
    std::vector<double> G{1., 0.1, 0.2, 0.5, 0.25};

    for (std::size_t i = 0; i < demog.get_set_deme_sizes().size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(demog.get_set_growth_rates()[i].when, when[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_growth_rates()[i].deme, deme[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_growth_rates()[i].G, G[i]);
        }
}

BOOST_AUTO_TEST_CASE(test_selfing_changes)
{
    std::vector<std::uint32_t> when{5, 5, 5, 10, 10};
    std::vector<int> deme{0, 6, 3, 1, 0};
    std::vector<double> S{1., 0., 0.2, 0.5, 0.25};

    for (std::size_t i = 0; i < demog.get_set_selfing_rates().size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(demog.get_set_selfing_rates()[i].when, when[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_selfing_rates()[i].deme, deme[i]);
            BOOST_REQUIRE_EQUAL(demog.get_set_selfing_rates()[i].S, S[i]);
        }
}

BOOST_AUTO_TEST_SUITE_END()
