#include <boost/test/unit_test.hpp>
#include <fwdpy11/regions/Sregion.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/policies/mutation.hpp>

BOOST_AUTO_TEST_SUITE(test_effect_sizes_of_zero)

// This test is motivated by GitHub issue #432.
// When a DES/DFE returns an effect size of zero,
// such variants were not put into the 'smutations'
// field of a genome, due to the way the fwdpp mutation
// base class constructor was being called.  In PR #433,
// we now require a DFE to tell the Mutation type if the
// variant is neutral or not, and the expectation is that
// the answer is "False".  This test creates a new Sregion
// type where all variants have effect size zero and we run
// a quick sim to make sure that they are all in both genomes
// and in the mutation table.

struct EsizeZero : public fwdpy11::Sregion
{
    EsizeZero(const fwdpy11::Region& r)
        : fwdpy11::Sregion(r, 1., 1, fwdpy11::process_input_dominance(1.0))
    {
    }

    std::unique_ptr<fwdpy11::Sregion>
    clone() const override
    {
        return std::unique_ptr<EsizeZero>(new EsizeZero(this->region));
    }

    double
    from_mvnorm(const double, const double) const override
    {
        return 0.0;
    }

    double
    generate_dominance(const fwdpy11::GSLrng_t&, const double) const override
    {
        return 1.;
    }

    std::uint32_t
    operator()(fwdpp::flagged_mutation_queue& recycling_bin,
               std::vector<fwdpy11::Mutation>& mutations,
               std::unordered_multimap<double, std::uint32_t>& lookup_table,
               const std::uint32_t generation,
               const fwdpy11::GSLrng_t& rng) const override
    {
        return fwdpy11::infsites_Mutation(
            recycling_bin, mutations, lookup_table, false, generation,
            [this, &rng]() { return region(rng); }, []() { return 0.; },
            [](const double) { return 1.; }, this->label());
    }
};

BOOST_AUTO_TEST_CASE(test_single_mutation)
{
    auto e = EsizeZero{fwdpy11::Region{0.0, 1.0, 1.0, false, 0}};
    auto pop = fwdpy11::DiploidPopulation(10, 1.0);
    auto q = fwdpp::flagged_mutation_queue{std::queue<std::size_t>{}};
    auto rng = fwdpy11::GSLrng_t(0);
    e(q, pop.mutations, pop.mut_lookup, pop.generation, rng);
    BOOST_REQUIRE_EQUAL(pop.mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(pop.mutations[0].neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
