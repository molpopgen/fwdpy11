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
    };
} // namespace fwdpy11

#endif
