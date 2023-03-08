#include <boost/test/unit_test.hpp>

#include <core/genetic_maps/regions.hpp>
#include <fwdpy11/rng.hpp>

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
