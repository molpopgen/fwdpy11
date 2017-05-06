#ifndef FWDPY11_UTIL_ADD_MUTATION_HPP__
#define FWDPY11_UTIL_ADD_MUTATION_HPP__

#include <fwdpy11/types.hpp>

std::size_t add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
                  const KTfwd::uint_t ncopies,
                  const std::tuple<double, double, double>& pos_s_h,
                  const std::uint16_t label);
std::size_t
add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::multilocus_t& pop,
             const std::size_t locus, const KTfwd::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label);
             
#endif
