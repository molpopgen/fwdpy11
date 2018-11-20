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
#ifndef FWDPY11_GENETIC_VALUES_DETAILS_GBR_HPP__
#define FWDPY11_GENETIC_VALUES_DETAILS_GBR_HPP__

#include <cmath>

namespace fwdpy11
{
    struct GBR
    /// "Gene-based recessice model of Thornton, Foran, Long 2013
    /// and Sanjak, Long, Thornton 2017
    {
        template <typename mutation_key_cont_t, typename mcont_t>
        inline double
        sum_haplotype_effect_sizes(const mutation_key_cont_t& keys,
                                   const mcont_t& mutations) const
        {
            double rv = 0.0;
            for (auto& k : keys)
                {
                    rv += mutations[k].s;
                }
            return rv;
        }
        template <typename diploid_t, typename gcont_t, typename mcont_t>
        inline double
        operator()(const diploid_t& dip, const gcont_t& gametes,
                   const mcont_t& mutations) const
        {
            double h1 = sum_haplotype_effect_sizes(
                gametes[dip.first].smutations, mutations);
            double h2 = sum_haplotype_effect_sizes(
                gametes[dip.second].smutations, mutations);
            return std::sqrt(h1 * h2);
        }
        GBR(){}
        GBR(const GBR &) = default;
    };
} // namespace fwdpy11

#endif
