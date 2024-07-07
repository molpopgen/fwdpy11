#include "fwdpp/fundamental_types/haploid_genome.hpp"
#include "fwdpp/fundamental_types/mutation_base.hpp"
#include "fwdpy11/genetic_value_data/genetic_value_data.hpp"
#include "fwdpy11/genetic_values/site_dependent_genetic_value.hpp"
#include "fwdpy11/types/Diploid.hpp"
#include "fwdpy11/types/DiploidPopulation.hpp"
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <fwdpy11/genetic_values/DiploidMultiplicative.hpp>
#include <fwdpy11/genetic_values/DiploidAdditive.hpp>

struct Mutation : public fwdpp::mutation_base
{
    double s;
    double h;
    Mutation(double pos, double s, double h)
        : fwdpp::mutation_base(pos, 0, false), s(s), h(h)
    {
    }
};

struct Pop
{
    std::vector<Mutation> mutations;
    std::vector<fwdpp::haploid_genome> genomes;
    std::vector<fwdpy11::DiploidMetadata> metadata;

    Pop() : mutations{}, genomes{}, metadata{}
    {
    }
};

struct Fwdpy11Pop
{
    fwdpy11::DiploidPopulation pop;
    std::vector<fwdpy11::DiploidMetadata> offspring_metadata;

    Fwdpy11Pop() : pop{1, 100.}, offspring_metadata()
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_low_level_gvalue_calculation)

BOOST_FIXTURE_TEST_CASE(test_multiplicative_fitness_0, Pop)
{
    mutations.emplace_back(1.0, -0.55, 1.0);
    genomes.emplace_back(1);
    genomes.emplace_back(1);
    genomes[0].smutations.emplace_back(0);
    BOOST_REQUIRE_EQUAL(genomes.size(), 2);
    fwdpy11::site_dependent_genetic_value sdgv{[](const double) { return false; }};
    auto gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w *= (1. + 2.0 * m.s); },
        [](double &w, const auto &m) { w *= (1. + m.s * m.h); },
        [](const double w) { return std::max(w, 0.0); }, 1.0);
    BOOST_REQUIRE_EQUAL(gv, 1.0 - 0.55);
}

BOOST_FIXTURE_TEST_CASE(test_multiplicative_fitness_1, Pop)
{
    mutations.emplace_back(1.0, -1.1, 1.0);
    genomes.emplace_back(1);
    genomes.emplace_back(1);
    genomes[0].smutations.emplace_back(0);
    BOOST_REQUIRE_EQUAL(genomes.size(), 2);
    fwdpy11::site_dependent_genetic_value sdgv{[](const double) { return false; }};
    auto gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w *= (1. + 2.0 * m.s); },
        [](double &w, const auto &m) { w *= (1. + m.s * m.h); },
        [](const double w) { return std::max(w, 0.0); }, 1.0);
    BOOST_REQUIRE_EQUAL(gv, 0.0);
}

BOOST_FIXTURE_TEST_CASE(test_multiplicative_fitness_2, Pop)
{
    mutations.emplace_back(1.0, -1.1, 1.0);
    mutations.emplace_back(1.0, -1.1, 1.0);
    genomes.emplace_back(1);
    genomes.emplace_back(1);
    genomes[0].smutations.emplace_back(0);
    genomes[0].smutations.emplace_back(1);
    BOOST_REQUIRE_EQUAL(genomes.size(), 2);
    // Gives the previous behavior where two "more than lethal"
    // variants give w > 0.0 when multiplied together
    fwdpy11::site_dependent_genetic_value sdgv{[](const double) { return false; }};
    auto gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w *= (1. + 2.0 * m.s); },
        [](double &w, const auto &m) { w *= (1. + m.s * m.h); },
        fwdpy11::final_multiplicative_fitness{}, 1.0);
    BOOST_REQUIRE(gv > 0.0);

    // Redo with a new closure that will return a fitness of 0.0
    // when we first see w <= 0.0
    sdgv
        = fwdpy11::site_dependent_genetic_value{[](const double w) { return w <= 0.0; }};
    gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w *= (1. + 2.0 * m.s); },
        [](double &w, const auto &m) { w *= (1. + m.s * m.h); },
        fwdpy11::final_multiplicative_fitness{}, 1.0);
    BOOST_REQUIRE_EQUAL(gv, 0.0);
}

BOOST_FIXTURE_TEST_CASE(test_additive_fitness, Pop)
{
    mutations.emplace_back(1.0, -1.1, 1.0);
    mutations.emplace_back(1.0, -1.1, 1.0);
    genomes.emplace_back(1);
    genomes.emplace_back(1);
    genomes[0].smutations.emplace_back(0);
    genomes[0].smutations.emplace_back(1);
    BOOST_REQUIRE_EQUAL(genomes.size(), 2);
    // Gives the previous behavior where two "more than lethal"
    // variants give w > 0.0 when multiplied together
    fwdpy11::site_dependent_genetic_value sdgv{[](const double) { return false; }};
    auto gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w += (2.0 * m.s); },
        [](double &w, const auto &m) { w += (m.s * m.h); },
        fwdpy11::final_additive_fitness(), 0.0);
    BOOST_REQUIRE_EQUAL(gv, 0.0);

    // Additive models DIFFER!!
    // If we use the same callback as for the multiplicative case,
    // we incorrectly get a final fitness > 0.
    sdgv
        = fwdpy11::site_dependent_genetic_value{[](const double w) { return w <= 0.0; }};
    gv = sdgv(
        genomes[0].smutations.begin(), genomes[0].smutations.end(),
        genomes[1].smutations.begin(), genomes[1].smutations.end(), mutations,
        [](double &w, const auto &m) { w += (m.s); },
        [](double &w, const auto &m) { w += (m.s * m.h); },
        fwdpy11::final_additive_fitness(), 0.0);
    BOOST_REQUIRE(gv > 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

// NOTE: the tests below rely on knowledge of how
// code that is not really part of the public API
// works and are thus fragile.

BOOST_AUTO_TEST_SUITE(test_diploid_multiplicative)

BOOST_FIXTURE_TEST_CASE(test_multiplicative_fitness_0, Fwdpy11Pop)
{
    pop.mutations.emplace_back(false, 0.1, -0.55, 1.0, 0);
    pop.haploid_genomes[0].smutations.emplace_back(0);
    pop.haploid_genomes.emplace_back(1);
    BOOST_REQUIRE_EQUAL(pop.haploid_genomes.size(), 2);
    auto m = fwdpy11::multiplicative_fitness_model(1, 2, nullptr);
    offspring_metadata.emplace_back(
        fwdpy11::DiploidMetadata{0., 0., 1., {0., 0., 0.}, 0, {0, 0}, 0, 0, {0, 0}});
    fwdpy11::GSLrng_t rng(42);
    pop.diploids[0] = {0, 1};
    fwdpy11::DiploidGeneticValueData data{rng,
                                          pop,
                                          pop.diploid_metadata[0],
                                          pop.diploid_metadata[0],
                                          0,
                                          offspring_metadata[0]};
    m(data);
    BOOST_REQUIRE_CLOSE(offspring_metadata[0].g, 1.0 - 0.55, 1e-3);
    BOOST_REQUIRE_CLOSE(offspring_metadata[0].w, 1.0 - 0.55, 1e-3);
}

BOOST_FIXTURE_TEST_CASE(test_multiplicative_fitness_1, Fwdpy11Pop)
{
    pop.mutations.emplace_back(false, 0.1, -1.55, 1.0, 0);
    pop.mutations.emplace_back(false, 0.2, -1.55, 1.0, 0);
    pop.haploid_genomes[0].smutations.emplace_back(0);
    pop.haploid_genomes[0].smutations.emplace_back(1);
    pop.haploid_genomes.emplace_back(1);
    BOOST_REQUIRE_EQUAL(pop.haploid_genomes.size(), 2);
    // This lambda will NOT exit when w first is <= 0.0
    fwdpy11::DiploidMultiplicative m(
        1, 2., fwdpy11::final_multiplicative_fitness(),
        [](const double) { return false; }, nullptr, nullptr);
    offspring_metadata.emplace_back(
        fwdpy11::DiploidMetadata{0., 0., 1., {0., 0., 0.}, 0, {0, 0}, 0, 0, {0, 0}});
    fwdpy11::GSLrng_t rng(42);
    pop.diploids[0] = {0, 1};
    BOOST_CHECK_EQUAL(pop.haploid_genomes[0].smutations.size(), 2);
    fwdpy11::DiploidGeneticValueData data{rng,
                                          pop,
                                          pop.diploid_metadata[0],
                                          pop.diploid_metadata[0],
                                          0,
                                          offspring_metadata[0]};
    m(data);
    BOOST_REQUIRE(offspring_metadata[0].g > 0.0);
    BOOST_REQUIRE(offspring_metadata[0].w > 0.0);

    // we now pass in this correct lambda expression
    m = fwdpy11::multiplicative_fitness_model(1, 2., nullptr);
    m(data);
    BOOST_REQUIRE_EQUAL(offspring_metadata[0].g, 0.);
    BOOST_REQUIRE_EQUAL(offspring_metadata[0].w, 0.);
}

BOOST_AUTO_TEST_SUITE_END()
