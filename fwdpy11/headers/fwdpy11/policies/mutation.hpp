//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef FWDPY11_POLICIES_HPP__
#define FWDPY11_POLICIES_HPP__

#include <cstdint>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpp/simfunctions/recycling.hpp>

namespace fwdpy11
{
    template <typename position_function, typename effect_size_function,
              typename dominance_function>
    std::size_t
    infsites_Mutation(fwdpp::flagged_mutation_queue &recycling_bin,
                      Population::mcont_t &mutations,
                      Population::lookup_table_t &lookup,
                      const fwdpp::uint_t &generation,
                      const position_function &posmaker,
                      const effect_size_function &esize_maker,
                      const dominance_function &hmaker,
                      const decltype(Mutation::xtra) x = 0)
    /*!
     * Mutation function to add a fwdpp::Mutation to a population.
     *
     * Generate mutations with single effect size.
     *
     * This implementation is a modification of fwdpp::infsites_popgenmut,
     * from the fwdpp library
     *
     * In order to use this function, it must be bound to a callable
     * that is a valid mutation function.  See examples for details.
     *
     * \param recycling_bin Recycling queue for mutations.
     * \param mutations Container of mutations
     * \param lookup Lookup table for mutation positions
     * \param generation The generation that is being mutated
     * \param pselected  The probability that a new mutation affects fitness
     * \param posmaker A function generating a mutation position.  Must be
     * convertible to std::function<double()>.
     * \param esize_maker A function to generate an effect size, given that a
     * mutation affects fitness. Must be convertible to
     * std::function<double()>.
     * \param hmaker A function to generate a dominance value, given that a
     * mutation affects fitness. Must be convertible to
     * std::function<double()>.
     *
     * \note "Neutral" mutations get assigned a dominance of zero.  The xtra
     * field is not written to.
     *
     */
    {
        auto pos = posmaker();
        while (lookup.find(pos) != lookup.end())
            {
                pos = posmaker();
            }
        auto idx = fwdpp::recycle_mutation_helper(recycling_bin, mutations,
                                                  pos, esize_maker(), hmaker(),
                                                  generation, x);
        lookup.emplace(pos, idx);
        return idx;
    }

    template <typename position_function, typename fixed_effect_size_function,
              typename fixed_dominance_function,
              typename effect_sizes_function,
              typename dominance_values_function>
    std::size_t
    infsites_Mutation(fwdpp::flagged_mutation_queue &recycling_bin,
                      Population::mcont_t &mutations,
                      Population::lookup_table_t &lookup,
                      const fwdpp::uint_t &generation,
                      const position_function &posmaker,
                      const fixed_effect_size_function &fixed_esize_maker,
                      const fixed_dominance_function &fixed_hmaker,
                      const effect_sizes_function &esizes,
                      const dominance_values_function &dominance,
                      const decltype(Mutation::xtra) x = 0)
    /*!
     * Mutation function to add a fwdpp::Mutation to a population.
     *
     * Generate mutations with single effect size.
     *
     * This implementation is a modification of fwdpp::infsites_popgenmut,
     * from the fwdpp library
     *
     * In order to use this function, it must be bound to a callable
     * that is a valid mutation function.  See examples for details.
     *
     * \param recycling_bin Recycling queue for mutations.
     * \param mutations Container of mutations
     * \param lookup Lookup table for mutation positions
     * \param generation The generation that is being mutated
     * \param pselected  The probability that a new mutation affects fitness
     * \param posmaker A function generating a mutation position.  Must be
     * convertible to std::function<double()>.
     * \param esize_maker A function to generate an effect size, given that a
     * mutation affects fitness. Must be convertible to
     * std::function<double()>.
     * \param hmaker A function to generate a dominance value, given that a
     * mutation affects fitness. Must be convertible to
     * std::function<double()>.
     *
     * \note "Neutral" mutations get assigned a dominance of zero.  The xtra
     * field is not written to.
     *
     */
    {
        auto pos = posmaker();
        while (lookup.find(pos) != lookup.end())
            {
                pos = posmaker();
            }
        auto idx = fwdpp::recycle_mutation_helper(
            recycling_bin, mutations, pos, fixed_esize_maker(), fixed_hmaker(),
            generation, esizes(), dominance(), x);
        lookup.emplace(pos, idx);
        return idx;
    }
} // namespace fwdpy11

#endif
