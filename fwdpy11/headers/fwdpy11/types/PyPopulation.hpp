#ifndef FWDPY11_PYPOPULATION_HPP__
#define FWDPY11_PYPOPULATION_HPP__

#include <tuple>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include <fwdpp/sugar/matrix.hpp>
#include "../rng.hpp"
#include "Diploid.hpp"

namespace fwdpy11
{
    template <typename mutation_type, typename mcont, typename gcont,
              typename mvector, typename ftvector, typename lookup_table_type>
    class PyPopulation
        : public fwdpp::sugar::popbase<mutation_type, mcont, gcont, mvector,
                                       ftvector, lookup_table_type>
    // Abstract base class (ABC) for population types
    {
      private:
        virtual void process_individual_input() = 0;

      public:
        using fwdpp_base
            = fwdpp::sugar::popbase<mutation_type, mcont, gcont, mvector,
                                    ftvector, lookup_table_type>;
        fwdpp::uint_t N;
        fwdpp::uint_t generation;

        std::vector<DiploidMetadata> diploid_metadata;

        virtual ~PyPopulation() = default;

        PyPopulation(PyPopulation &&) = default;
        PyPopulation(const PyPopulation &) = default;

        PyPopulation(fwdpp::uint_t N_)
            : fwdpp_base{ N_ }, N{ N_ }, generation{ 0 }, diploid_metadata(N)
        {
        }

        template <typename gametes_input, typename mutations_input>
        explicit PyPopulation(
            const fwdpp::uint_t N_, gametes_input &&g, mutations_input &&m,
            typename fwdpp_base::gamete_t::mutation_container::size_type
                reserve_size)
            : fwdpp_base{ std::forward<gametes_input>(g),
                          std::forward<mutations_input>(m), reserve_size },
              N{ N_ }, generation{ 0 }, diploid_metadata(N)
        {
        }

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &gametes)
            = 0;

        std::int64_t
        find_mutation_by_key(
            const std::tuple<double, double, fwdpp::uint_t> &key,
            const std::int64_t offset) const
        {
            auto itr = std::find_if(
                this->mutations.begin() + offset, this->mutations.end(),
                [&key](const typename mvector::value_type &mutation) {
                    return key
                           == std::tie(mutation.pos, mutation.s, mutation.g);
                });
            if (itr == this->mutations.end())
                return -1;
            return static_cast<std::int64_t>(
                std::distance(this->mutations.begin(), itr));
        }

        std::int64_t
        find_fixation_by_key(
            const std::tuple<double, double, fwdpp::uint_t> &key,
            const std::int64_t offset) const
        {
            auto itr = std::find_if(
                this->fixations.begin() + offset, this->fixations.end(),
                [&key](const typename mvector::value_type &mutation) {
                    return key
                           == std::tie(mutation.pos, mutation.s, mutation.g);
                });
            if (itr == this->fixations.end())
                return -1;
            return static_cast<std::int64_t>(
                std::distance(this->fixations.begin(), itr));
        }

        virtual fwdpp::data_matrix
        sample_individuals(const std::vector<std::size_t> &,
                           const bool) const = 0;
        virtual fwdpp::data_matrix
        sample_random_individuals(const GSLrng_t &, const std::uint32_t,
                                  const bool) const = 0;

        // The next two functions can be called from derived
        // classes to help implementing the above two functions.
        template <typename poptype>
        fwdpp::data_matrix
        sample_individuals_details(const poptype &pop,
                                   const std::vector<std::size_t> &individuals,
                                   const bool haplotype) const
        {
            if (std::any_of(individuals.begin(), individuals.end(),
                            [&pop](const std::size_t i) {
                                return i >= pop.diploids.size();
                            }))
                {
                    throw std::out_of_range("individual index out of range");
                }
            auto keys = fwdpp::mutation_keys(pop, individuals, true, true);
            using vtype = decltype(keys.first);
            using key_type = typename vtype::value_type;
            const auto sorting_lambda = [&pop](const key_type &a,
                                               const key_type &b) {
                return pop.mutations[a.first].pos < pop.mutations[b.first].pos;
            };
            //sort keys on position
            std::sort(keys.first.begin(), keys.first.end(), sorting_lambda);
            std::sort(keys.second.begin(), keys.second.end(), sorting_lambda);

            //Remove keys to mutations that are fixed in the sample
            keys.first.erase(
                std::remove_if(keys.first.begin(), keys.first.end(),
                               [&individuals](const key_type &k) {
                                   return k.second == 2 * individuals.size();
                               }),
                keys.first.end());
            keys.second.erase(
                std::remove_if(keys.second.begin(), keys.second.end(),
                               [&individuals](const key_type &k) {
                                   return k.second == 2 * individuals.size();
                               }),
                keys.second.end());
            return (haplotype == true)
                       ? fwdpp::haplotype_matrix(pop, individuals, keys.first,
                                                 keys.second)
                       : fwdpp::genotype_matrix(pop, individuals, keys.first,
                                                keys.second);
        }

        template <typename poptype>
        fwdpp::data_matrix
        sample_random_individuals_details(const poptype &pop,
                                          const GSLrng_t &rng,
                                          const std::uint32_t nsam,
                                          const bool haplotype) const
        {
            if (nsam > pop.N)
                {
                    throw std::invalid_argument(
                        "sample size > population size");
                }
            std::vector<std::size_t> all_individuals(pop.N);
            std::iota(all_individuals.begin(), all_individuals.end(), 0);
            if (nsam == pop.N)
                {
                    return sample_individuals_details(pop, all_individuals,
                                                      haplotype);
                }

            std::vector<std::size_t> individuals(nsam, 0);
            gsl_ran_choose(rng.get(), individuals.data(), individuals.size(),
                           all_individuals.data(), all_individuals.size(),
                           sizeof(std::size_t));
            return sample_individuals_details(pop, individuals, haplotype);
        }
    };
} // namespace fwdpy11

#endif
