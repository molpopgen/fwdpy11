#pragma once

#include <fwdpp/type_traits.hpp>
#include <limits>

namespace fwdpy11
{
    struct site_dependent_genetic_value
    {
        std::function<bool(const double)> clamp;
        //! The return value type
        using result_type = double;

        template <typename iterator_t, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het,
                  typename make_return_value>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Range-based call operator.  Calculates genetic values over ranges of
          mutation keys first1/last1
          and first2/last2, which are iterators derived from the 'smutations'
          of two haploid_genomes in a diploid.

          \param first1 Iterator to first mutation derived from haploid_genome 1
          \param last1 Iterator to one past the last mutation derived from
          haploid_genome 1
          \param first2 Iterator to first mutation derived from haploid_genome 2
          \param last2 Iterator to one past the last mutation derived from
          haploid_genome 2
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param make_return_value Policy generated the final return value.
                 Must be equivalent to std::function<double(double)>
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            static_assert(fwdpp::traits::is_mutation<
                              typename MutationContainerType::value_type>::value,
                          "MutationContainerType::value_type must be a mutation type");
            static_assert(
                std::is_convertible<
                    updating_policy_hom,
                    std::function<void(double &, const typename MutationContainerType::
                                                     value_type &)>>::value,
                "decltype(fpol_hom) must be convertible to "
                "std::function<void(double &,const typename "
                "MutationContainerType::value_type");
            static_assert(
                std::is_convertible<
                    updating_policy_het,
                    std::function<void(double &, const typename MutationContainerType::
                                                     value_type &)>>::value,
                "decltype(fpol_het) must be convertible to "
                "std::function<void(double &,const typename "
                "MutationContainerType::value_type");
            result_type w = starting_value;
            if (first1 == last1 && first2 == last2)
                return rv_function(w);
            else if (first1 == last1)
                {
                    for (; first2 != last2; ++first2)
                        {
                            fpol_het(w, mutations[*first2]);
                            if (clamp(w))
                                {
                                    return rv_function(0.0);
                                }
                        }
                    return rv_function(w);
                }
            else if (first2 == last2)
                {
                    for (; first1 != last1; ++first1)
                        {
                            fpol_het(w, mutations[*first1]);
                            if (clamp(w))
                                {
                                    return rv_function(0.0);
                                }
                        }

                    return rv_function(w);
                }
            for (; first1 != last1; ++first1)
                {
                    for (; first2 != last2 && *first1 != *first2
                           && mutations[*first2].pos < mutations[*first1].pos;
                         ++first2)
                        // All mutations in this range are Aa
                        {
                            fpol_het(w, mutations[*first2]);
                            if (clamp(w))
                                {
                                    return rv_function(0.0);
                                }
                        }
                    if (first2 < last2
                        && (*first1 == *first2
                            || mutations[*first1].pos == mutations[*first2].pos))
                        // mutation with index first1 is homozygous
                        {
                            fpol_hom(w, mutations[*first1]);
                            if (clamp(w))
                                {
                                    return rv_function(0.0);
                                }
                            ++first2; // increment so that we don't re-process
                            // this site as a het next time 'round
                        }
                    else // mutation first1 is heterozygous
                        {
                            fpol_het(w, mutations[*first1]);
                            if (clamp(w))
                                {
                                    return rv_function(0.0);
                                }
                        }
                }
            for (; first2 != last2; ++first2)
                {
                    fpol_het(w, mutations[*first2]);
                    if (clamp(w))
                        {
                            return rv_function(0.0);
                        }
                }
            return rv_function(w);
        }

        template <typename iterator_t, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                first1, last1, first2, last2, mutations, fpol_hom, fpol_het,
                [](double d) { return d; }, starting_value);
        }

        template <typename HaploidGenomeType, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het,
                  typename make_return_value>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Calculates genetic value for a diploid whose genotype
          across sites is given by haploid_genomes g1 and  g2.

          \param g1 A haploid_genome
          \param g2 A haploid_genome
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            static_assert(fwdpp::traits::is_haploid_genome<HaploidGenomeType>::value,
                          "HaploidGenomeType::value_type must be a haploid_genome "
                          "type");
            return this->operator()(g1.smutations.cbegin(), g1.smutations.cend(),
                                    g2.smutations.cbegin(), g2.smutations.cend(),
                                    mutations, fpol_hom, fpol_het, rv_function,
                                    starting_value);
        }

        template <typename HaploidGenomeType, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                g1, g2, mutations, fpol_hom, fpol_het, [](double d) { return d; },
                starting_value);
        }

        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType, typename updating_policy_hom,
                  typename updating_policy_het, typename make_return_value>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Calculates genetic value for a diploid type.

          \param dip A diploid
          \param haploid_genomes The container of haploid_genomes for the simulation.
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            return this->operator()(haploid_genomes[dip.first],
                                    haploid_genomes[dip.second], mutations, fpol_hom,
                                    fpol_het, rv_function, starting_value);
        }

        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType, typename updating_policy_hom,
                  typename updating_policy_het>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                dip, haploid_genomes, mutations, fpol_hom, fpol_het,
                [](double d) { return d; }, starting_value);
        }
    };
}
