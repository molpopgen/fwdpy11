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
#include <fwdpp/ts/std_table_collection.hpp>
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
        static std::shared_ptr<fwdpp::ts::std_table_collection>
        init_tables(const fwdpp::uint_t N, const double L)
        // Default initialization of tables.
        // We ensure that there are 2N nodes in pop zero
        // at time 0
        {
            if (L == std::numeric_limits<double>::max())
                {
                    return std::make_shared<fwdpp::ts::std_table_collection>(L);
                }
            return std::make_shared<fwdpp::ts::std_table_collection>(2 * N, 0, 0, L);
            // return fwdpp::ts::std_table_collection(2 * N, 0, 0, L);
        }

        bool
        tables_equal(const Population &rhs) const
        {
            // This is correct/validated as of fwdpp 0.7.2
            return *tables == *rhs.tables;
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
        bool is_simulating;

        std::shared_ptr<fwdpp::ts::std_table_collection> tables;
        std::vector<fwdpp::ts::table_index_t> alive_nodes, preserved_sample_nodes;

        // These track genetic values for the individuals differently
        // from what diploid_metadata/ancient_sample_metadata do.
        // For the case of a multivariate trait, we imagine this to
        // represent a matrix of N rows by "dimensions" columns.
        std::vector<double> genetic_value_matrix, ancient_sample_genetic_value_matrix;

        Population(fwdpp::uint_t N_, const double L)
            : fwdpp_base{N_}, N{N_}, generation{0}, is_simulating{false},
              tables(init_tables(N_, L)), alive_nodes{}, preserved_sample_nodes{},
              genetic_value_matrix{}, ancient_sample_genetic_value_matrix{}
        {
        }

        virtual ~Population() = default;
        Population(Population &&) = default;
        Population(const Population &) = default;
        Population &operator=(const Population &) = default;
        Population &operator=(Population &&) = default;

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

        void
        rebuild_mutation_lookup(bool from_tables)
        {
            this->mut_lookup.clear();
            if (from_tables)
                {
                    for (const auto &mr : this->tables->mutations)
                        {
                            if (mr.key >= this->mutations.size())
                                {
                                    throw std::runtime_error(
                                        "rebuild_mutation_lookup: mutation "
                                        "table key out of range");
                                }
                            this->mut_lookup.emplace(
                                this->tables->sites[mr.site].position, mr.key);
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
        virtual void fill_preserved_nodes() = 0;
    };
}

#endif
