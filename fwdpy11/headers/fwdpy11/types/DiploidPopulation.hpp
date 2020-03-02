#ifndef FWDPY11_POP_HPP__
#define FWDPY11_POP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include "create_pops.hpp"
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <fwdpp/poptypes/tags.hpp>
#include <fwdpp/sugar/add_mutation.hpp>

namespace fwdpy11
{
    class DiploidPopulation : public Population
    {
      private:
        void
        process_individual_input()
        {
            std::vector<fwdpp::uint_t> gcounts(this->haploid_genomes.size(),
                                               0);
            for (auto &&dip : diploids)
                {
                    this->validate_individual_keys(dip.first);
                    this->validate_individual_keys(dip.second);
                    gcounts[dip.first]++;
                    gcounts[dip.second]++;
                }
            this->validate_haploid_genome_counts(gcounts);
        }

        void
        fill_sample_nodes_from_metadata(
            std::vector<fwdpp::ts::TS_NODE_INT> &n,
            const std::vector<fwdpy11::DiploidMetadata> &metadata)
        {
            n.clear();
            for (auto &md : metadata)
                {
                    n.push_back(md.nodes[0]);
                    n.push_back(md.nodes[1]);
                }
        }

        void
        finish_construction(const std::vector<std::uint32_t> &deme_sizes)
        {
            if (!this->N)
                {
                    throw std::invalid_argument("population size must be > 0");
                }
            std::size_t label = 0;
            for (auto &&d : this->diploid_metadata)
                {
                    d.label = label++;
                    d.w = 1.0;
                    if (this->tables.genome_length()
                        == std::numeric_limits<double>::max())
                        {
                            d.nodes[0] = d.nodes[1] = -1;
                        }
                    else // Fix GitHub issue 332
                        {
                            d.nodes[0] = 2 * d.label;
                            d.nodes[1] = 2 * d.label + 1;
                        }
                }
            if (this->tables.genome_length()
                != std::numeric_limits<double>::max())
                {
                    using itype = decltype(this->diploid_metadata[0].nodes[0]);
                    using ntype = std::remove_reference<itype>::type;
                    if (label >= std::numeric_limits<ntype>::max())
                        {
                            throw std::invalid_argument(
                                "population size too large for node "
                                "type");
                        }
                    std::int32_t deme_label = 0;
                    std::size_t i = 0;
                    for (auto ds : deme_sizes)
                        {
                            for (std::uint32_t j = 0; j < ds; ++j, ++i)
                                {
                                    auto &md = diploid_metadata[i];
                                    md.deme = deme_label;
                                    tables.node_table[md.nodes[0]].deme
                                        = deme_label;
                                    tables.node_table[md.nodes[1]].deme
                                        = deme_label;
                                }
                            ++deme_label;
                        }
                }
        }

      public:
        using dipvector_t = std::vector<DiploidGenotype>;
        using diploid_t = dipvector_t::value_type;
        using popbase_t = Population;
        using popmodel_t = fwdpp::poptypes::DIPLOID_TAG;
        using fitness_t
            = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                           typename popbase_t::gcont_t,
                                           typename popbase_t::mcont_t>;

        dipvector_t diploids;
        //TODO: initalized ancient_sample_metadata and tables,
        //TODO figure out what to do with class constructor??
        //TODO Introduce types for ancient sample individual and node tracking
        std::vector<DiploidMetadata> diploid_metadata, ancient_sample_metadata;
        std::vector<ancient_sample_record> ancient_sample_records;

        // Constructors for Python
        DiploidPopulation(const fwdpp::uint_t N, const double length)
            : Population{ N, length }, diploids(N, { 0, 0 }),
              diploid_metadata(N), ancient_sample_metadata{},
              ancient_sample_records{}
        {
            finish_construction({ N });
        }

        DiploidPopulation(const std::vector<std::uint32_t> &deme_sizes,
                          const double length)
            : Population{ std::accumulate(begin(deme_sizes), end(deme_sizes),
                                          0u),
                          length },
              diploids(std::accumulate(begin(deme_sizes), end(deme_sizes), 0u),
                       { 0, 0 }),
              diploid_metadata(N), ancient_sample_metadata{},
              ancient_sample_records{}

        {
            finish_construction(deme_sizes);
        }

        template <typename diploids_input, typename genomes_input,
                  typename mutations_input>
        explicit DiploidPopulation(diploids_input &&d, genomes_input &&g,
                                   mutations_input &&m)
            : Population(static_cast<fwdpp::uint_t>(d.size()),
                         std::forward<genomes_input>(g),
                         std::forward<mutations_input>(m), 100),
              diploids(std::forward<diploids_input>(d)), diploid_metadata(N),
              ancient_sample_metadata{}, ancient_sample_records{}
        //! Constructor for pre-determined population status
        {
            this->process_individual_input();
        }

        ~DiploidPopulation() = default;
        DiploidPopulation(DiploidPopulation &&) = default;
        DiploidPopulation(const DiploidPopulation &) = default;
        DiploidPopulation &operator=(const DiploidPopulation &) = default;
        DiploidPopulation &operator=(DiploidPopulation &&) = default;

        virtual bool
        operator==(const DiploidPopulation &rhs) const
        {
            return this->diploids == rhs.diploids
                   && this->diploid_metadata == rhs.diploid_metadata
                   && this->ancient_sample_metadata
                          == rhs.ancient_sample_metadata
                   && popbase_t::test_equality(rhs);
        };

        void
        clear()
        {
            diploids.clear();
            popbase_t::clear_containers();
        }

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &new_mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &haploid_genomes)
        {
            std::unordered_set<double> poschecker;
            for (const auto &m : new_mutations)
                {
                    if (this->mut_lookup.find(m.pos) != this->mut_lookup.end())
                        {
                            throw std::invalid_argument(
                                "attempting to add new mutation at "
                                "already-mutated position");
                        }
                    if (poschecker.find(m.pos) != poschecker.end())
                        {
                            throw std::invalid_argument(
                                "attempting to add multiple mutations at the "
                                "same position");
                        }
                    poschecker.insert(m.pos);
                }
            std::vector<std::size_t> rv;

            for (auto &i : new_mutations)
                {
                    auto pos = i.pos;
                    // remaining preconditions get checked by fwdpp:
                    auto idx = fwdpp::add_mutation(
                        *this, individuals, haploid_genomes, std::move(i));

                    // fwdpp's function doesn't update the lookup:
                    this->mut_lookup.emplace(pos, idx);
                    rv.push_back(idx);
                }
            return rv;
        }

        fwdpp::data_matrix
        sample_individuals(const std::vector<std::size_t> &individuals,
                           const bool haplotype, const bool remove_fixed) const
        {
            return sample_individuals_details(*this, individuals, haplotype,
                                              remove_fixed);
        }

        fwdpp::data_matrix
        sample_random_individuals(const GSLrng_t &rng,
                                  const std::uint32_t nsam,
                                  const bool haplotype,
                                  const bool remove_fixed) const
        {
            return sample_random_individuals_details(*this, rng, nsam,
                                                     haplotype, remove_fixed);
        }

        void
        fill_alive_nodes() override
        {
            fill_sample_nodes_from_metadata(alive_nodes, diploid_metadata);
        }

        void
        fill_preserved_nodes()
        {
            fill_sample_nodes_from_metadata(preserved_sample_nodes,
                                            ancient_sample_metadata);
        }

        virtual std::size_t
        ancient_sample_metadata_size() const override
        {
            return ancient_sample_metadata.size();
        }
    };
} // namespace fwdpy11
#endif
