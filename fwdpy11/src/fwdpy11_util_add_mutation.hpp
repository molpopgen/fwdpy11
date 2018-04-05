#ifndef FWDPY11_UTIL_ADD_MUTATION_HPP__
#define FWDPY11_UTIL_ADD_MUTATION_HPP__

#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

std::size_t add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::SlocusPop& pop,
                  const fwdpp::uint_t ncopies,
                  const std::tuple<double, double, double>& pos_s_h,
                  const std::uint16_t label);
std::size_t
add_mutation(const fwdpy11::GSLrng_t& rng, fwdpy11::MlocusPop& pop,
             const std::size_t locus, const fwdpp::uint_t ncopies,
             const std::tuple<double, double, double>& pos_s_h,
             const std::uint16_t label);
             
#endif
