#ifndef FWDPY11_PYPOPULATION_HPP__
#define FWDPY11_PYPOPULATION_HPP__

#include <tuple>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/data_matrix.hpp>
#include <fwdpp/ts/table_collection.hpp>
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

      public:
        using fwdpp_base
            = fwdpp::sugar::popbase<mutation_type, mcont, gcont, mvector,
                                    ftvector, lookup_table_type>;
        fwdpp::uint_t N;
        fwdpp::uint_t generation;

        //TODO: initalized ancient_sample_metadata and tables,
        //TODO figure out what to do with class constructor??
        //TODO Introduce types for ancient sample individual and node tracking
        std::vector<DiploidMetadata> diploid_metadata, ancient_sample_metadata;
        std::vector<ancient_sample_record> ancient_sample_records;
        std::vector<fwdpp::uint_t> mcounts_from_preserved_nodes;
        fwdpp::ts::table_collection tables;

        virtual ~PyPopulation() = default;

        PyPopulation(PyPopulation &&) = default;
        PyPopulation(const PyPopulation &) = default;

        PyPopulation(fwdpp::uint_t N_, const double L)
            : fwdpp_base{ N_ }, N{ N_ }, generation{ 0 },
              diploid_metadata(N), ancient_sample_metadata{},
              ancient_sample_records{}, mcounts_from_preserved_nodes{},
              tables(init_tables(N_, L))
        {
        }

        template <typename gametes_input, typename mutations_input>
        explicit PyPopulation(
            const fwdpp::uint_t N_, gametes_input &&g, mutations_input &&m,
            typename fwdpp_base::gamete_t::mutation_container::size_type
                reserve_size)
            : fwdpp_base{ std::forward<gametes_input>(g),
                          std::forward<mutations_input>(m), reserve_size },
              N{ N_ }, generation{ 0 },
              diploid_metadata(N), ancient_sample_metadata{},
              ancient_sample_records{}, mcounts_from_preserved_nodes{},
              tables(std::numeric_limits<double>::max())
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
        sample_individuals(const std::vector<std::size_t> &, const bool,
                           const bool) const = 0;
        virtual fwdpp::data_matrix
        sample_random_individuals(const GSLrng_t &, const std::uint32_t,
                                  const bool, const bool) const = 0;

        // The next two functions can be called from derived
        // classes to help implementing the above two functions.
        template <typename poptype>
        fwdpp::data_matrix
        sample_individuals_details(const poptype &pop,
                                   const std::vector<std::size_t> &individuals,
                                   const bool haplotype,
                                   const bool remove_fixed) const
        {
            if (std::any_of(individuals.begin(), individuals.end(),
                            [&pop](const std::size_t i) {
                                return i >= pop.diploids.size();
                            }))
                {
                    throw std::out_of_range("individual index out of range");
                }
            if (haplotype)
                {
                    // returns a haplotype matrix
                    return fwdpp::sample_individuals(pop, individuals, true,
                                                     true, remove_fixed);
                }
            auto keys = fwdpp::fwdpp_internal::generate_filter_sort_keys(
                pop, individuals, true, true, remove_fixed);
            return fwdpp::genotype_matrix(pop, individuals, keys.first,
                                          keys.second);
        }

        template <typename poptype>
        fwdpp::data_matrix
        sample_random_individuals_details(const poptype &pop,
                                          const GSLrng_t &rng,
                                          const std::uint32_t nsam,
                                          const bool haplotype,
                                          const bool remove_fixed) const
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
                                                      haplotype, remove_fixed);
                }

            std::vector<std::size_t> individuals(nsam, 0);
            gsl_ran_choose(rng.get(), individuals.data(), individuals.size(),
                           all_individuals.data(), all_individuals.size(),
                           sizeof(std::size_t));
            return sample_individuals_details(pop, individuals, haplotype,
                                              remove_fixed);
        }

        // TODO: redo once fwdpp adds comparison operator for table_collection
        bool
        tables_equal(const PyPopulation &rhs) const
        {
            return mcounts_from_preserved_nodes
                   == rhs.mcounts_from_preserved_nodes;
            if (tables.genome_length() != rhs.tables.genome_length())
                {
                    return false;
                }
            if (tables.edge_table.size() != rhs.tables.edge_table.size()
                || tables.node_table.size() != rhs.tables.node_table.size()
                || tables.mutation_table.size()
                       != rhs.tables.mutation_table.size())
                {
                    return false;
                }
            for (std::size_t i = 0; i < tables.edge_table.size(); ++i)
                {
                    auto e1 = tables.edge_table[i];
                    auto e2 = rhs.tables.edge_table[i];
                    bool l = (e1.left == e2.left);
                    bool r = (e1.right == e2.right);
                    bool p = (e1.parent == e2.parent);
                    bool c = (e1.child == e2.child);
                    if (!l || !r || !p || !c)
                        {
                            return false;
                        }
                }
            for (std::size_t i = 0; i < tables.node_table.size(); ++i)
                {
                    auto n1 = tables.node_table[i];
                    auto n2 = rhs.tables.node_table[i];
                    bool p = (n1.population == n2.population);
                    bool t = (n1.time == n2.time);
                    if (!p || !t)
                        {
                            return false;
                        }
                }
            for (std::size_t i = 0; i < tables.mutation_table.size(); ++i)
                {
                    auto mr1 = tables.mutation_table[i];
                    auto mr2 = rhs.tables.mutation_table[i];
                    bool n = (mr1.node == mr2.node);
                    bool k = (mr1.key == mr2.key);
                    if (!n || !k)
                        {
                            return false;
                        }
                }
        }
    };
} // namespace fwdpy11

#endif
