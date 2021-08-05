#include <boost/test/unit_test.hpp>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <fwdpy11/discrete_demography/MigrationMatrix.hpp>
#include <fwdpy11/discrete_demography/constants.hpp>
#include <fwdpy11/discrete_demography/simulation/mating_event_type.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include "discrete_demography_roundtrips.hpp"
#include "discrete_demography_util.hpp"
#include "discrete_demography_fixtures.hpp"

using namespace fwdpy11::discrete_demography;

using deme_sizes_t = std::unordered_map<int, int>;

BOOST_FIXTURE_TEST_SUITE(test_lowlevel_demographic_events, population_fixture)

BOOST_AUTO_TEST_CASE(test_simple_moves_from_single_deme)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 1);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    BOOST_CHECK_EQUAL(deme_sizes.size(), 2);
    for (auto&& i : deme_sizes)
        {
            BOOST_CHECK_EQUAL(i.second, 50);
        }
}

BOOST_AUTO_TEST_CASE(test_simple_copies_from_single_deme)
{
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 0.5, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 1);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    BOOST_CHECK_EQUAL(deme_sizes.size(), 2);
    BOOST_CHECK_EQUAL(deme_sizes[0], 100);
    BOOST_CHECK_EQUAL(deme_sizes[1], 50);
}

BOOST_AUTO_TEST_CASE(test_simple_back_and_forth_copy)
{
    mass_migrations.emplace_back(copy_individuals(1, 0, 1, 0.5, true));
    mass_migrations.emplace_back(copy_individuals(2, 1, 0, 1.0, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 3);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    BOOST_CHECK_EQUAL(deme_sizes.size(), 2);
    deme_sizes_t expected{{0, 150}, {1, 50}};
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_simple_copies_from_multiple_demes)
{
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    mass_migrations.emplace_back(copy_individuals(1, 1, 2, 0.5, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 2);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    deme_sizes_t expected{{0, 50}, {1, 50}, {2, 25}};
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_single_deme_growth)
{
    unsigned N1 = 3412;
    unsigned t = 111;
    double G = std::exp((std::log(N1) - std::log(pop.N)) / t);
    set_growth_rates.emplace_back(16, 0, G);
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 15 + t + 1);
    BOOST_CHECK_EQUAL(pop.N, N1);
    BOOST_CHECK_EQUAL(pop.diploid_metadata.size(), static_cast<std::size_t>(N1));
}

BOOST_AUTO_TEST_CASE(test_two_deme_growth)
{
    std::vector<unsigned> N0{90, 10};
    std::vector<unsigned> t{14, 23};
    std::vector<unsigned> N1{5361, 616};
    auto G0 = std::exp((std::log(N1[0]) - std::log(N0[0])) / t[0]);
    auto G1 = std::exp((std::log(N1[1]) - std::log(N0[1])) / t[1]);

    set_growth_rates.emplace_back(7, 0, G0);
    set_growth_rates.emplace_back(7 + t[0], 0, fwdpy11::discrete_demography::NOGROWTH);
    set_growth_rates.emplace_back(33, 1, G1);
    set_growth_rates.emplace_back(33 + t[1], 1, fwdpy11::discrete_demography::NOGROWTH);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.1, true));

    auto demog = make_model();

    DiscreteDemography_roundtrip(rng, pop, demog, 100);
    BOOST_CHECK_EQUAL(pop.N, N1[0] + N1[1]);
}

BOOST_AUTO_TEST_CASE(test_two_deme_growth_with_hard_reset)
{
    std::vector<unsigned> N0{90, 10};
    std::vector<unsigned> t{14, 23};
    std::vector<unsigned> N1{5361, 616};
    auto G0 = std::exp((std::log(N1[0]) - std::log(N0[0])) / t[0]);
    auto G1 = std::exp((std::log(N1[1]) - std::log(N0[1])) / t[1]);

    set_growth_rates.emplace_back(7, 0, G0);
    set_growth_rates.emplace_back(33, 1, G1);
    set_growth_rates.emplace_back(33 + t[1], 1, fwdpy11::discrete_demography::NOGROWTH);

    // Cut off the growth in deme 0 after a few generations,
    // and manually set the new deme size to 100 w/no growth
    set_deme_sizes.emplace_back(11, 0, 100, true);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.1, true));

    auto demog = make_model();

    DiscreteDemography_roundtrip(rng, pop, demog, 100);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    decltype(N1) expected = {100, N1[1]};
    BOOST_CHECK_EQUAL(pop.N, expected[0] + expected[1]);
    for (auto&& ds : deme_sizes)
        {
            BOOST_CHECK_EQUAL(ds.second, expected[static_cast<std::size_t>(ds.first)]);
        }
}

BOOST_AUTO_TEST_CASE(test_two_deme_growth_without_hard_reset)
{
    std::vector<unsigned> N0{90, 10};
    std::vector<unsigned> t{14, 23};
    std::vector<unsigned> N1{5361, 616};
    auto G0 = std::exp((std::log(N1[0]) - std::log(N0[0])) / t[0]);
    auto G1 = std::exp((std::log(N1[1]) - std::log(N0[1])) / t[1]);

    set_growth_rates.emplace_back(7, 0, G0);
    set_growth_rates.emplace_back(7 + t[0], 0, fwdpy11::discrete_demography::NOGROWTH);
    set_growth_rates.emplace_back(33, 1, G1);
    set_growth_rates.emplace_back(33 + t[1], 1, fwdpy11::discrete_demography::NOGROWTH);

    // Cut off the growth in deme 0 after a few generations,
    // and manually set the new deme size to 100 AND CONTINUE GROWTH
    set_deme_sizes.emplace_back(11, 0, 100, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.1, true));

    auto demog = make_model();

    DiscreteDemography_roundtrip(rng, pop, demog, 100);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    decltype(N1) expected = {
        static_cast<unsigned>(std::round(100.0 * std::pow(G0, 7 + t[0] - 11))), N1[1]};
    BOOST_CHECK_EQUAL(pop.N, expected[0] + expected[1]);
    for (auto&& ds : deme_sizes)
        {
            BOOST_CHECK_EQUAL(ds.second, expected[static_cast<std::size_t>(ds.first)]);
        }
}

BOOST_AUTO_TEST_CASE(test_two_moves_in_same_generation)
/*
 * In generation 0, deme 0 splits equally
 * into demes 1 and 2.  This should leave
 * deme 0 empty and the sizes of demes
 * 1 and 2 both equal to 0.5 the initial
 * size.
 */
{
    deme_sizes_t expected = {{1, pop.N / 2}, {2, pop.N / 2}};
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 1);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_copies_happen_before_moves)
/*
 * This is a more complex test.
 * In the same generation:
 * 1. All individuals are copied from deme 0 to deme 1.
 * 2. 50% of deme 0 moves to deme 2.
 * 3. 50% of deme 0 moves to deme 3.
 * 
 * This, at the end, the size of deme 1 should be
 * the initial size and the size of demes 2 and 3
 * should be 0.5*initial size
 */
{
    deme_sizes_t expected = {{1, pop.N}, {2, pop.N / 2}, {3, pop.N / 2}};
    mass_migrations.emplace_back(move_individuals(0, 0, 3, 0.5, true));
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 1, true));
    mass_migrations.emplace_back(move_individuals(0, 0, 2, 0.5, true));
    auto demog = make_model();
    DiscreteDemography_roundtrip(rng, pop, demog, 1);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_mass_move_with_growth)
/*
 * In generation 5, the mass movement from 0 to 1
 * will reset the growth rate in deme 0 to
 * fwdpy11.NOGROWTH
 * 
 * This test is also handy b/c growth results in
 * an odd total N by the time the mass migration happens
 */
{
    mass_migrations.emplace_back(move_individuals(5, 0, 1, 0.5, true));
    double G = 1.1;
    set_growth_rates.emplace_back(0, 0, G);
    auto ddemog = make_model();
    auto N0{pop.N};
    DiscreteDemography_roundtrip(rng, pop, ddemog, 10);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    auto N5 = std::round(N0 * std::pow(G, 5));
    // NOTE: the cast to int gives us Python's // here, kinda.
    deme_sizes_t expected{{0, N5 / 2}, {1, N5 - static_cast<int>(N5 / 2)}};
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_mass_move_with_growth_no_reset)
{
    mass_migrations.emplace_back(move_individuals(5, 0, 1, 0.5, false));
    double G = 1.1;
    set_growth_rates.emplace_back(0, 0, G);
    auto ddemog = make_model();
    auto N0{pop.N};
    DiscreteDemography_roundtrip(rng, pop, ddemog, 10);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    auto N5 = std::round(N0 * std::pow(G, 5));
    auto N_after_mass_mig_0 = std::round(static_cast<int>(N5 / 2) * std::pow(G, 5));
    deme_sizes_t expected{{0, N_after_mass_mig_0}, {1, N5 - static_cast<int>(N5 / 2)}};
    for (auto&& e : expected)
        {
            BOOST_CHECK(deme_sizes.find(e.first) != end(deme_sizes));
            BOOST_CHECK_EQUAL(e.second, deme_sizes[e.first]);
        }
}

BOOST_AUTO_TEST_CASE(test_simple_two_deme_migration)
{
    set_migmatrix(std::vector<double>{0.5, 0.5, 0.5, 0.5}, 2, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    auto ddemog = make_model();
    auto N{pop.N};
    DiscreteDemography_roundtrip(rng, pop, ddemog, 5);
    BOOST_CHECK_EQUAL(N, pop.N);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    BOOST_CHECK_EQUAL(deme_sizes.size(), 2);
    for (auto&& ds : deme_sizes)
        {
            BOOST_CHECK_EQUAL(ds.second, pop.N / 2);
        }
}

BOOST_AUTO_TEST_CASE(test_change_migration_rates_simple_two_deme_migration)
/*
 *  For a 2-deme model, the mig matrix is
 * [[0, 1]
 *  [1, 0]]
 * so that all offspring have both parents from the other deme,
 * which gives us an easy check on how many migration events
 * will be recorded by the test simulation.
 * 
 * After 3 generations, we reset the migration rates to be
 * [[1, 0],
 *  [1, 0]]
 * so that all parents are from deme zero.
 */
{
    set_migmatrix(std::vector<double>{0., 1., 1., 0}, 2, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    set_migration_rates.emplace_back(3, std::vector<double>{1, 0, 1, 0});
    auto demog = make_model();
    auto N{pop.N};
    auto migevents = DiscreteDemography_roundtrip(rng, pop, demog, 5);
    BOOST_CHECK_EQUAL(N, pop.N);
    BOOST_CHECK_EQUAL(migevents.size(), 2 * N * 3 + 2 * N);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    for (auto&& ds : deme_sizes)
        {
            BOOST_CHECK_EQUAL(ds.second, N / 2);
        }
}

BOOST_AUTO_TEST_CASE(test_selfing_vs_migration)
/*
 * Parents of deme 0 are all migrants from deme 1.
 * Parents of deme 1 are all migrants from deme 0.
 * Deme 0 never selfs.  Deme 1 always selfs.
 */
{
    set_migmatrix(std::vector<double>{0., 1., 1., 0}, 2, false);
    mass_migrations.emplace_back(move_individuals(0, 0, 1, 0.5, true));
    set_selfing_rates.emplace_back(0, 1, 1.);
    auto demog = make_model();
    auto migevents = DiscreteDemography_roundtrip(rng, pop, demog, 5);
    for (auto& m : migevents)
        {
            if (m.parental_deme == 1) // parent from deme 1
                {
                    BOOST_CHECK(
                        m.mating_event
                        == fwdpy11::discrete_demography::mating_event_type::selfing);
                    BOOST_CHECK(m.parental_deme != m.offspring_deme);
                }
            else if (m.parental_deme == 0)
                {
                    BOOST_CHECK(
                        m.mating_event
                        == fwdpy11::discrete_demography::mating_event_type::outcrossing);
                    BOOST_CHECK(m.parental_deme != m.offspring_deme);
                }
            else
                {
                    throw std::runtime_error("unexpected result");
                }
        }
}

BOOST_AUTO_TEST_CASE(test_ghost_population)
{
    mass_migrations.emplace_back(copy_individuals(5, 0, 1, 1., true));
    mass_migrations.emplace_back(copy_individuals(10, 1, 0, 0.1, true));
    set_deme_sizes.emplace_back(10, 1, 0, true);
    set_deme_sizes.emplace_back(10, 0, 100, true);
    auto ddemog = make_model();
    DiscreteDemography_roundtrip(rng, pop, ddemog, 20);
    BOOST_CHECK_EQUAL(pop.N, 100);
    auto deme_sizes = get_deme_sizes(pop.diploid_metadata);
    BOOST_CHECK_EQUAL(deme_sizes.size(), 1);
    BOOST_CHECK(deme_sizes.find(0) != end(deme_sizes));
    BOOST_CHECK(deme_sizes.find(1) == end(deme_sizes));
    BOOST_CHECK_EQUAL(deme_sizes[0], 100);
}

BOOST_AUTO_TEST_CASE(test_exponential_decline)
{
    set_growth_rates.emplace_back(0, 0, 0.5);
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      fwdpy11::discrete_demography::GlobalExtinction);
}

// NOTE: the tests below are pretty extreme.
// They are working hard to trigger exceptions

BOOST_AUTO_TEST_CASE(bad_metadata_label_when_mass_migration_happens)
{
    mass_migrations.emplace_back(copy_individuals(0, 0, 1, 1., true));
    // individual's array index no longer matches label. BAD
    pop.diploid_metadata[0].label += 7;
    auto ddemog = make_model();
    BOOST_CHECK_THROW({ DiscreteDemography_roundtrip(rng, pop, ddemog, 20); },
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(mass_migration_via_change_in_migration_rates)
{
    // Set deme size 0 to 0 immediately.
    set_deme_sizes.emplace_back(0, 0, 0, true);
    // Create a new deme, that will get 100% of ancestry from
    // individuals in deme 0
    set_deme_sizes.emplace_back(0, 1, pop.N, true);

    // Set an initial migration matrix representing just 100%
    // ancestry of ancestral pop to itself.
    set_migmatrix(std::vector<double>{1., 0., 0., 0.}, 2, false);

    // Change the migration matrix to reflect the "mass migration"
    // event founding deme 1
    set_migration_rates.emplace_back(0, std::vector<double>{0., 0., 1., 0.});

    // The mass migration event is done, so all the ancestry
    // must now come from the new deme.
    set_migration_rates.emplace_back(1, std::vector<double>{0., 0., 0., 1.});

    auto demog = make_model();
    auto migevents = DiscreteDemography_roundtrip(rng, pop, demog, 5);

    // We only record migrations if parent deme != offspring deme,
    // so all migrations should have generation 1 (the first offspring
    // generation) and be from 0 -> 1, representing the pulse migration.
    for (auto& m : migevents)
        {
            BOOST_REQUIRE_EQUAL(m.generation, 1);
            BOOST_REQUIRE_EQUAL(m.parental_deme, 0);
            BOOST_REQUIRE_EQUAL(m.offspring_deme, 1);
        }
}

BOOST_AUTO_TEST_SUITE_END()

