#ifndef FWDPY11_POPULATION_HPP__
#define FWDPY11_POPULATION_HPP__

#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include <unordered_map>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/poptypes/popbase.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/data_matrix.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include "../rng.hpp"
#include "Mutation.hpp"

namespace fwdpy11
{
    class Population
        : public fwdpp::poptypes::popbase<
              Mutation, std::vector<Mutation>, std::vector<fwdpp::haploid_genome>,
              std::vector<Mutation>, std::vector<fwdpp::uint_t>,
              std::unordered_multimap<double, fwdpp::uint_t>>
    // Base class for population types
    {
      private:
        virtual void process_individual_input() = 0;

        static fwdpp::ts::table_collection
        init_tables(const fwdpp::uint_t N, const double L)
        // Default initialization of tables.
        // We ensure that there are 2N nodes in pop zero
        // at time 0
        {
            if (L == std::numeric_limits<double>::max())
                {
                    return fwdpp::ts::table_collection(L);
                }
            return fwdpp::ts::table_collection(2 * N, 0, 0, L);
        }

        bool
        tables_equal(const Population &rhs) const
        {
            // This is correct/validated as of fwdpp 0.7.2
            return tables == rhs.tables;
        }

      public:
        using mutation_type = Mutation;
        using mutation_vector = std::vector<mutation_type>;
        using genome_type = fwdpp::haploid_genome;
        using genome_vector = std::vector<genome_type>;
        using mutation_position_hash = std::unordered_multimap<double, fwdpp::uint_t>;
        using fwdpp_base
            = fwdpp::poptypes::popbase<mutation_type, mutation_vector, genome_vector,
                                       mutation_vector, std::vector<fwdpp::uint_t>,
                                       mutation_position_hash>;

        fwdpp::uint_t N;
        fwdpp::uint_t generation;

        fwdpp::ts::table_collection tables;
        std::vector<fwdpp::ts::TS_NODE_INT> alive_nodes, preserved_sample_nodes;

        // These track genetic values for the individuals differently
        // from what diploid_metadata/ancient_sample_metadata do.
        // For the case of a multivariate trait, we imagine this to
        // represent a matrix of N rows by "dimensions" columns.
        std::vector<double> genetic_value_matrix, ancient_sample_genetic_value_matrix;

        Population(fwdpp::uint_t N_, const double L)
            : fwdpp_base{N_}, N{N_}, generation{0},
              tables(init_tables(N_, L)), alive_nodes{}, preserved_sample_nodes{},
              genetic_value_matrix{}, ancient_sample_genetic_value_matrix{}
        {
        }

        template <typename gametes_input, typename mutations_input>
        explicit Population(
            const fwdpp::uint_t N_, gametes_input &&g, mutations_input &&m,
            typename fwdpp_base::haploid_genome_t::mutation_container::size_type
                reserve_size)
            : fwdpp_base{std::forward<gametes_input>(g),
                         std::forward<mutations_input>(m), reserve_size},
              N{N_}, generation{0}, tables(std::numeric_limits<double>::max()),
              alive_nodes{}, preserved_sample_nodes{}, genetic_value_matrix{},
              ancient_sample_genetic_value_matrix{}
        {
        }

        virtual ~Population() = default;
        Population(Population &&) = default;
        Population(const Population &) = default;
        Population &operator=(const Population &) = default;
        Population &operator=(Population &&) = default;

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &gametes)
            = 0;

        std::int64_t
        find_mutation_by_key(const std::tuple<double, double, fwdpp::uint_t> &key,
                             const std::int64_t offset) const
        {
            auto itr = std::find_if(
                this->mutations.begin() + offset, this->mutations.end(),
                [&key](const typename mutation_vector::value_type &mutation) {
                    return key == std::tie(mutation.pos, mutation.s, mutation.g);
                });
            if (itr == this->mutations.end())
                return -1;
            return static_cast<std::int64_t>(
                std::distance(this->mutations.begin(), itr));
        }

        std::int64_t
        find_fixation_by_key(const std::tuple<double, double, fwdpp::uint_t> &key,
                             const std::int64_t offset) const
        {
            auto itr = std::find_if(
                this->fixations.begin() + offset, this->fixations.end(),
                [&key](const typename mutation_vector::value_type &mutation) {
                    return key == std::tie(mutation.pos, mutation.s, mutation.g);
                });
            if (itr == this->fixations.end())
                return -1;
            return static_cast<std::int64_t>(
                std::distance(this->fixations.begin(), itr));
        }

        virtual fwdpp::data_matrix sample_individuals(const std::vector<std::size_t> &,
                                                      const bool, const bool) const = 0;
        virtual fwdpp::data_matrix sample_random_individuals(const GSLrng_t &,
                                                             const std::uint32_t,
                                                             const bool,
                                                             const bool) const = 0;

        // The next two functions can be called from derived
        // classes to help implementing the above two functions.
        template <typename poptype>
        fwdpp::data_matrix
        sample_individuals_details(const poptype &pop,
                                   const std::vector<std::size_t> &individuals,
                                   const bool haplotype, const bool remove_fixed) const
        {
            if (std::any_of(
                    individuals.begin(), individuals.end(),
                    [&pop](const std::size_t i) { return i >= pop.diploids.size(); }))
                {
                    throw std::out_of_range("individual index out of range");
                }
            if (haplotype)
                {
                    // returns a haplotype matrix
                    return fwdpp::sample_individuals(pop, individuals, true, true,
                                                     remove_fixed);
                }
            auto keys = fwdpp::fwdpp_internal::generate_filter_sort_keys(
                pop, individuals, true, true, remove_fixed);
            return fwdpp::genotype_matrix(pop, individuals, keys.first, keys.second);
        }

        template <typename poptype>
        fwdpp::data_matrix
        sample_random_individuals_details(const poptype &pop, const GSLrng_t &rng,
                                          const std::uint32_t nsam, const bool haplotype,
                                          const bool remove_fixed) const
        {
            if (nsam > pop.N)
                {
                    throw std::invalid_argument("sample size > population size");
                }
            std::vector<std::size_t> all_individuals(pop.N);
            std::iota(all_individuals.begin(), all_individuals.end(), 0);
            if (nsam == pop.N)
                {
                    return sample_individuals_details(pop, all_individuals, haplotype,
                                                      remove_fixed);
                }

            std::vector<std::size_t> individuals(nsam, 0);
            gsl_ran_choose(rng.get(), individuals.data(), individuals.size(),
                           all_individuals.data(), all_individuals.size(),
                           sizeof(std::size_t));
            return sample_individuals_details(pop, individuals, haplotype, remove_fixed);
        }

        void
        rebuild_mutation_lookup(bool from_tables)
        {
            this->mut_lookup.clear();
            if (from_tables)
                {
                    for (const auto &mr : this->tables.mutation_table)
                        {
                            if (mr.key >= this->mutations.size())
                                {
                                    throw std::runtime_error(
                                        "rebuild_mutation_lookup: mutation "
                                        "table key out of range");
                                }
                            this->mut_lookup.emplace(
                                this->tables.site_table[mr.site].position, mr.key);
                        }
                }
            else
                {
                    if (this->mutations.size() != this->mcounts.size()
                        || (!this->mcounts_from_preserved_nodes.empty()
                            && (this->mutations.size()
                                != this->mcounts_from_preserved_nodes.size())))
                        {
                            throw std::runtime_error(
                                "rebuild_mutation_lookup: container size "
                                "mismatch");
                        }
                    for (std::size_t i = 0; i < this->mutations.size(); ++i)
                        {
                            if (this->mcounts[i] > 0
                                || (!this->mcounts_from_preserved_nodes.empty()
                                    && this->mcounts_from_preserved_nodes[i] > 0))
                                {
                                    this->mut_lookup.emplace(this->mutations[i].pos, i);
                                }
                        }
                }
        }

        bool
        test_equality(const Population &rhs) const
        {
            return this->is_equal(rhs) && tables_equal(rhs)
                   && this->genetic_value_matrix == rhs.genetic_value_matrix
                   && this->ancient_sample_genetic_value_matrix
                          == rhs.ancient_sample_genetic_value_matrix;
        }

        virtual std::size_t ancient_sample_metadata_size() const = 0;
        virtual void fill_alive_nodes() = 0;
    };
}

#endif
