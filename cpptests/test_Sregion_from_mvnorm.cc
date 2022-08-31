#include <boost/test/unit_test.hpp>
#include <cmath>
#include <fwdpy11/regions/ExpS.hpp>
#include <fwdpy11/regions/GaussianS.hpp>
#include <fwdpy11/regions/UniformS.hpp>
#include <fwdpy11/regions/GammaS.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

static const auto P = 0.2;

auto
make_region()
{
    return fwdpy11::Region(0., 1., 1.0, false, 0);
}

BOOST_AUTO_TEST_SUITE(test_Sregion_from_mvnorm)

BOOST_AUTO_TEST_CASE(test_ExpS)
{
    auto mean = 0.2;
    auto e = fwdpy11::ExpS(make_region(), 1.0, mean, fwdpy11::fixed_dominance(1.));
    auto d = e.from_mvnorm(0, P);
    BOOST_CHECK(std::isfinite(d));
    BOOST_CHECK_CLOSE(P, gsl_cdf_exponential_P(d, mean), 1e-8);
}

BOOST_AUTO_TEST_CASE(test_GaussianS)
{
    auto sd = 0.2;
    auto e = fwdpy11::GaussianS(make_region(), 1.0, sd, fwdpy11::fixed_dominance(1.));
    auto d = e.from_mvnorm(0, P);
    BOOST_CHECK(std::isfinite(d));
    BOOST_CHECK_CLOSE(P, gsl_cdf_gaussian_P(d, sd), 1e-8);
}

BOOST_AUTO_TEST_CASE(test_GammaS)
{
    auto a = 0.2;
    auto b = 10;
    double mean = b * a;
    auto e = fwdpy11::GammaS(make_region(), 1.0, mean, a, fwdpy11::fixed_dominance(1.));
    auto d = e.from_mvnorm(0, P);
    BOOST_CHECK(std::isfinite(d));
    BOOST_CHECK_CLOSE(P, gsl_cdf_gamma_P(d, a, b), 1e-8);
}

BOOST_AUTO_TEST_CASE(test_UniformS)
{
    auto lo = 0.1;
    auto hi = 0.47;
    auto e = fwdpy11::UniformS(make_region(), 1.0, lo, hi, fwdpy11::fixed_dominance(1.));
    auto d = e.from_mvnorm(0, P);
    BOOST_CHECK(std::isfinite(d));
    BOOST_CHECK_CLOSE(P, gsl_cdf_flat_P(d, lo, hi), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
