#include <cstdint>
#include <queue>
#include <vector>
#include <unordered_map>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_matrix.h>
#include <fwdpy11/regions/ExpS.hpp>
#include <fwdpy11/regions/mvDES.hpp>
#include <fwdpy11/regions/LogNormalS.hpp>
#include <fwdpy11/regions/MultivariateGaussianEffects.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpy11/regions/Region.hpp>
#include <fwdpy11/regions/Sregion.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

BOOST_AUTO_TEST_SUITE(test_MutationDominance)

static fwdpy11::Region
make_dummy_region()
{
    return fwdpy11::Region(0, 1, 1, true, 0);
}

static auto
generate_mutation(const fwdpy11::Sregion &s, std::vector<fwdpy11::Mutation> &mutations)
{
    fwdpy11::GSLrng_t rng(42);
    fwdpp::flagged_mutation_queue q(std::queue<std::size_t>{});
    std::unordered_multimap<double, std::uint32_t> lookup_table;
    return s(q, mutations, lookup_table, 0, rng);
}

BOOST_AUTO_TEST_CASE(test_fixed_dominance_init_and_clone)
{
    fwdpy11::GSLrng_t rng(42);
    auto f = fwdpy11::fixed_dominance(1.0);
    BOOST_CHECK_EQUAL(f(rng, 0.2135123), 1.);

    auto fc(f);
    BOOST_CHECK_EQUAL(fc(rng, 0.2135123), 1.);
}

BOOST_AUTO_TEST_CASE(test_fixed_dominance_round_trip)
{
    fwdpy11::ExpS e(make_dummy_region(), 2, 1., fwdpy11::fixed_dominance(0.25));
    std::vector<fwdpy11::Mutation> mutations;
    auto x = generate_mutation(e, mutations);
    BOOST_CHECK_EQUAL(x, 0);
    BOOST_CHECK_EQUAL(mutations[x].h, 0.25);
}

BOOST_AUTO_TEST_CASE(test_fixed_dominance_round_trip_mvDES_vector_of_Sregion)
{
    std::vector<std::unique_ptr<fwdpy11::Sregion>> odist;
    for (int i = 0; i < 3; ++i)
        {
            odist.emplace_back(std::unique_ptr<fwdpy11::Sregion>(new fwdpy11::ExpS(
                make_dummy_region(), 2, 1., fwdpy11::fixed_dominance(i))));
        }
    std::vector<double> vcov(9, 0.);
    gsl_matrix_view vcov_view = gsl_matrix_view_array(vcov.data(), 3, 3);
    gsl_matrix_set_identity(&vcov_view.matrix);
    fwdpy11::mvDES mv(odist, std::vector<double>(3, 0.), vcov_view.matrix);
    std::vector<fwdpy11::Mutation> mutations;
    auto x = generate_mutation(mv, mutations);
    BOOST_CHECK_EQUAL(x, 0);
    for (int i = 0; i < 3; ++i)
        {
            BOOST_CHECK_EQUAL(mutations[x].heffects[i], static_cast<double>(i));
        }
}

BOOST_AUTO_TEST_CASE(test_fixed_dominance_round_trip_mvDES_LogNormalS)
{
    fwdpy11::LogNormalS l(make_dummy_region(), 1, fwdpy11::fixed_dominance(1. / 8.));
    std::vector<double> vcov(9, 0.);
    gsl_matrix_view vcov_view = gsl_matrix_view_array(vcov.data(), 3, 3);
    gsl_matrix_set_identity(&vcov_view.matrix);
    fwdpy11::mvDES mv(l, std::vector<double>(3, 0.), vcov_view.matrix);
    std::vector<fwdpy11::Mutation> mutations;
    auto x = generate_mutation(mv, mutations);
    BOOST_CHECK_EQUAL(x, 0);
    for (int i = 0; i < 3; ++i)
        {
            BOOST_CHECK_EQUAL(mutations[x].heffects[i], 1. / 8.);
        }
}

BOOST_AUTO_TEST_CASE(test_fixed_dominance_round_trip_mvDES_MultivariateGaussian)
{
    std::vector<double> vcov(9, 0.);
    gsl_matrix_view vcov_view = gsl_matrix_view_array(vcov.data(), 3, 3);
    gsl_matrix_set_identity(&vcov_view.matrix);
    fwdpy11::MultivariateGaussianEffects mvg(make_dummy_region(), 2., vcov_view.matrix,
                                             1, fwdpy11::fixed_dominance(1. / 6.));
    fwdpy11::mvDES mv(mvg, std::vector<double>(3, 0.));
    std::vector<fwdpy11::Mutation> mutations;
    auto x = generate_mutation(mv, mutations);
    BOOST_CHECK_EQUAL(x, 0);
    for (int i = 0; i < 3; ++i)
        {
            BOOST_CHECK_EQUAL(mutations[x].heffects[i], 1. / 6.);
        }
}

BOOST_AUTO_TEST_SUITE_END()
