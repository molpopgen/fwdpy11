#include "fwdpy11_util_add_mutation.hpp"
#include <fwdpp/sugar/add_mutation.hpp>
#include <tuple>
#include <cstdint>
#include <array>
#include <numeric>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>

void
add_next_dip(const fwdpy11::GSLrng_t& rng, std::vector<std::size_t>& dips,
             const KTfwd::uint_t N)
{
    unsigned next_dip = gsl_rng_uniform_int(rng.get(), N);
    while (std::find(dips.begin(), dips.end(), next_dip) != dips.end())
        {
            next_dip = gsl_rng_uniform_int(rng.get(), N);
        }
    dips.push_back(next_dip);
}

std::tuple<std::vector<std::size_t>, std::vector<short>>
get_diploids_and_genos(const fwdpy11::GSLrng_t& rng, KTfwd::uint_t ncopies,
                       const KTfwd::uint_t N)
// Figure out who the diploids are that get the new mutation,
// and their corresponding genotype
{
    std::vector<std::size_t> dips;
    std::vector<short> genos;
    double q = static_cast<double>(ncopies) / static_cast<double>(2 * N);
    double pAa = 2. * (1. - q) * q;
    double paa = gsl_sf_pow_int(q, 2);
    pAa /= (pAa + paa); // cond'l prob of being Aa
    // Get # aa
    unsigned naa = gsl_ran_binomial(rng.get(), 1. - pAa, ncopies / 2);
    unsigned nAa = ncopies - 2 * naa;
    for (unsigned i = 0; i < naa; ++i)
        {
            genos.push_back(2);
            add_next_dip(rng, dips, N);
        }
    for (unsigned i = 0; i < nAa; ++i)
        {
            genos.push_back(1);
            add_next_dip(rng, dips, N);
        }
    if (std::accumulate(genos.begin(), genos.end(), 0u) != ncopies)
        {
            throw std::runtime_error(
                "fatal error: incorrect number of mutation copies generated");
        }
    return std::make_tuple(std::move(dips), std::move(genos));
}

std::size_t
add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
             const KTfwd::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label)
/// Add a new mutation at HW proportions
{
    if (pop.mut_lookup.find(std::get<0>(pos_s_h)) != pop.mut_lookup.end())
        {
            throw std::invalid_argument(
                "new mutation position already exists in population");
        }
    auto dips_and_genos = get_diploids_and_genos(rng, ncopies, pop.N);
    auto rv = KTfwd::add_mutation(pop, std::get<0>(dips_and_genos),
                                  std::get<1>(dips_and_genos),
                                  std::get<0>(pos_s_h), std::get<1>(pos_s_h),
                                  std::get<2>(pos_s_h), pop.generation, label);
    pop.mut_lookup.insert(pop.mutations[rv].pos);
    return rv;
}

std::size_t
add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::multilocus_t& pop,
             const std::size_t locus, const KTfwd::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label)
/// Add a new mutation at HW proportions
{
    if (pop.mut_lookup.find(std::get<0>(pos_s_h)) != pop.mut_lookup.end())
        {
            throw std::invalid_argument(
                "new mutation position already exists in population");
        }
    auto dips_and_genos = get_diploids_and_genos(rng, ncopies, pop.N);
    auto rv = KTfwd::add_mutation(pop, locus, std::get<0>(dips_and_genos),
                                  std::get<1>(dips_and_genos),
                                  std::get<0>(pos_s_h), std::get<1>(pos_s_h),
                                  std::get<2>(pos_s_h), pop.generation, label);
    pop.mut_lookup.insert(pop.mutations[rv].pos);
    return rv;
}
