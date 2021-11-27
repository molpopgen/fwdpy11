#include <vector>
#include <boost/test/unit_test.hpp>

#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GSSmo.hpp>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSSmo.hpp>

BOOST_AUTO_TEST_SUITE(test_gss_classes)

BOOST_AUTO_TEST_CASE(test_GSSmo_update)
{
    std::vector<fwdpy11::Optimum> optima{
        {0, 0.25, 1.0}, {10, 0.5, 1.0}, {11, 0.75, 1.0}, {12, 1.0, 1.0}, {50, 1.25, 1.0},
    };
    for (std::uint32_t g = 0; g < 100; ++g)
        {
            fwdpy11::DiploidPopulation pop(10000, 1.);
            pop.generation = g;
            fwdpy11::GSSmo gssmo(optima);
            for (std::uint32_t i = 0; i <= 51; ++i, ++pop.generation)
                {
                    gssmo.update(pop);
                    std::size_t j = 0;
                    while (j < optima.size() && optima[j].when <= pop.generation)
                        {
                            ++j;
                        }
                    BOOST_REQUIRE(j > 0);
                    BOOST_REQUIRE(j - 1 < optima.size());
                    BOOST_REQUIRE_EQUAL(optima[j - 1].opt, gssmo.opt);
                }
        }
}

BOOST_AUTO_TEST_CASE(test_MultivariateGSSmo_update)
{
    std::vector<fwdpy11::PleiotropicOptima> optima{
        {0, {0.25, 0.25}, 1.0}, {10, {0.5, 0.5}, 1.0},   {11, {0.75, 0.75}, 1.0},
        {12, {1.0, 1.0}, 1.0},  {50, {1.25, 1.25}, 1.0},
    };
    fwdpy11::DiploidPopulation pop(10000, 1.);
    fwdpy11::MultivariateGSSmo mvgssmo(optima);
    for (std::uint32_t i = 0; i <= 51; ++i, ++pop.generation)
        {
            mvgssmo.update(pop);
            std::size_t j = 0;
            while (j < optima.size() && optima[j].when <= pop.generation)
                {
                    ++j;
                }
            BOOST_REQUIRE(j > 0);
            BOOST_REQUIRE(j - 1 < optima.size());
            BOOST_REQUIRE(optima[j - 1].optima == mvgssmo.current_timepoint_optima);
        }
}

BOOST_AUTO_TEST_SUITE_END()
