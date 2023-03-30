#include "fwdpy11/regions/RecombinationRegions.hpp"
#include "fwdpy11/regions/Region.hpp"
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/unit_test_suite.hpp>
#include <core/genetic_maps/regions.hpp>
#include <fwdpy11/rng.hpp>
#include <limits>
#include <memory>

BOOST_AUTO_TEST_SUITE(test_core_genetic_map_regions)

BOOST_AUTO_TEST_CASE(test_binomial_region)
{
    fwdpy11::GSLrng_t rng(42);

    {
        fwdpy11_core::BinomialInterval b(0, 1, 0.5, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes > 0);
    }

    {
        fwdpy11_core::BinomialInterval b(0, 1, 0.0, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes == 0);
    }

    {
        fwdpy11_core::BinomialInterval b(0, 1, 1, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes == 100);
    }
}

BOOST_AUTO_TEST_CASE(test_binomial_point)
{
    fwdpy11::GSLrng_t rng(42);

    {
        fwdpy11_core::BinomialPoint b(0, 0.5, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes > 0);
    }

    {
        fwdpy11_core::BinomialPoint b(0, 0.0, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes == 0);
    }

    {
        fwdpy11_core::BinomialPoint b(0, 1, false);
        int successes = 0;
        for (int i = 0; i < 100; ++i)
            {
                std::vector<double> bp;
                b.breakpoint(rng, bp);
                if (!bp.empty())
                    {
                        successes++;
                    }
            }
        BOOST_CHECK(successes == 100);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_generalized_genetic_map)

BOOST_AUTO_TEST_CASE(test_binomial_point)
{
    std::vector<std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>> callbacks;
    callbacks.push_back(std::make_unique<fwdpy11_core::BinomialPoint>(0, 1, true));
    callbacks.push_back(std::make_unique<fwdpy11_core::BinomialPoint>(1, 1, true));
    auto map = fwdpy11::GeneralizedGeneticMap({}, std::move(callbacks));
    auto rng = fwdpy11::GSLrng_t(42);
    auto bp = map(rng);
    BOOST_REQUIRE_EQUAL(bp.size(), 3);
    BOOST_REQUIRE_EQUAL(bp[0], 0.0);
    BOOST_REQUIRE_EQUAL(bp[1], 1.0);
    BOOST_REQUIRE_EQUAL(bp[2], std::numeric_limits<double>::max());
}

BOOST_AUTO_TEST_CASE(test_poisson_point)
{
    std::vector<std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>> callbacks;
    // Mathematical near-certainty that there are some breakpoints
    callbacks.push_back(
        std::make_unique<fwdpy11_core::PoissonInterval>(0, 1, 1000.0, false));
    callbacks.push_back(
        std::make_unique<fwdpy11_core::PoissonInterval>(1, 2, 0.0, false));
    auto map = fwdpy11::GeneralizedGeneticMap(std::move(callbacks), {});
    auto rng = fwdpy11::GSLrng_t(42);
    auto bp = map(rng);
    BOOST_REQUIRE(!bp.empty());
    BOOST_REQUIRE_EQUAL(
        std::count_if(begin(bp), end(bp) - 1, [](auto x) { return x >= 0 && x < 1.0; }),
        bp.size() - 1);
}

BOOST_AUTO_TEST_CASE(test_binomial_interval_map)
{
    std::vector<fwdpy11::Region> regions;
    regions.push_back(fwdpy11::Region(0, 10, 1e-3, false, 0));
    regions.push_back(fwdpy11::Region(10, 20, 5e-3, true, 0));
    auto b = fwdpy11_core::BinomialIntervalMap(0.5, true, regions);
    BOOST_CHECK_EQUAL(b.left(), 0.0);
    BOOST_CHECK_EQUAL(b.right(), 20.0);
    std::vector<std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>> callbacks;
    callbacks.push_back(
        std::make_unique<fwdpy11_core::BinomialIntervalMap>(std::move(b)));
    auto map = fwdpy11::GeneralizedGeneticMap({}, std::move(callbacks));
    auto rng = fwdpy11::GSLrng_t(42);
    for (int i = 0; i < 100; ++i)
        {
            auto bp = map(rng);
        }
}

BOOST_AUTO_TEST_SUITE_END()
