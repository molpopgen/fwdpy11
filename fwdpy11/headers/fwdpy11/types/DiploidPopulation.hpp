#ifndef FWDPY11_POP_HPP__
#define FWDPY11_POP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include "fwdpy11/types/Mutation.hpp"
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <fwdpp/poptypes/tags.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/site_visitor.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>

namespace fwdpy11
{
    class DiploidPopulation : public Population
    {
      private:
        void
        fill_sample_nodes_from_metadata(
            std::vector<fwdpp::ts::table_index_t> &n,
            const std::vector<fwdpy11::DiploidMetadata> &metadata)
        {
            std::unordered_set<fwdpp::ts::table_index_t> unodes;
            const auto update = [&unodes, &n](auto i) {
                if (unodes.find(i) == end(unodes))
                    {
                        n.push_back(i);
                        unodes.insert(i);
                    }
            };
            n.clear();
            for (auto &md : metadata)
                {
                    for (auto i : md.nodes)
                        {
                            update(i);
                        }
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
                    if (this->tables->genome_length()
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
            if (this->tables->genome_length() != std::numeric_limits<double>::max())
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
                                    tables->nodes[md.nodes[0]].deme = deme_label;
                                    tables->nodes[md.nodes[1]].deme = deme_label;
                                }
                            ++deme_label;
                        }
                }
        }

      public:
        using dipvector_t = std::vector<DiploidGenotype>;
        using diploid_t = dipvector_t::value_type;
        using diploid_type = dipvector_t::value_type;
        using popbase_t = Population;
        using popmodel_t = fwdpp::poptypes::DIPLOID_TAG;
        using fitness_t
            = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                           typename popbase_t::genome_container,
                                           typename popbase_t::mutation_container>;

        dipvector_t diploids;
        //TODO: initalized ancient_sample_metadata and tables,
        //TODO figure out what to do with class constructor??
        //TODO Introduce types for ancient sample individual and node tracking
        std::vector<DiploidMetadata> diploid_metadata, ancient_sample_metadata;

        // Constructors for Python
        DiploidPopulation(const fwdpp::uint_t N, const double length)
            : Population{2, N, length}, diploids(N, {0, 0}),
              diploid_metadata(N), ancient_sample_metadata{}
        {
            finish_construction({N});
        }

        DiploidPopulation(const std::vector<std::uint32_t> &deme_sizes,
                          const double length)
            : Population{2, std::accumulate(begin(deme_sizes), end(deme_sizes), 0u),
                         length},
              diploids(std::accumulate(begin(deme_sizes), end(deme_sizes), 0u), {0, 0}),
              diploid_metadata(N), ancient_sample_metadata{}
        {
            finish_construction(deme_sizes);
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
                   && this->ancient_sample_metadata == rhs.ancient_sample_metadata
                   && popbase_t::test_equality(rhs);
        };

        void
        clear()
        {
            diploids.clear();
            popbase_t::clear_containers();
        }

        void
        fill_alive_nodes() override
        {
            fill_sample_nodes_from_metadata(alive_nodes, diploid_metadata);
        }

        void
        fill_preserved_nodes() override
        {
            fill_sample_nodes_from_metadata(preserved_sample_nodes,
                                            ancient_sample_metadata);
        }

        virtual std::size_t
        ancient_sample_metadata_size() const override
        {
            return ancient_sample_metadata.size();
        }

        virtual void
        record_ancient_samples(const std::vector<std::uint32_t> &individuals) override
        {
            for (auto individual : individuals)
                {
                    if (individual >= this->N)
                        {
                            throw std::invalid_argument(
                                "ancient sample index greater than "
                                "current population size");
                        }
                    ancient_sample_metadata.push_back(diploid_metadata[individual]);
                }
        }

        virtual void
        update_ancient_sample_genetic_value_matrix(
            const std::vector<std::uint32_t> &alive_individuals,
            std::size_t total_dim) override
        {
            for (auto individual : alive_individuals)
                {
                    if (individual >= this->N)
                        {
                            throw std::invalid_argument(
                                "ancient sample index greater than "
                                "current population size");
                        }
                    // Update the genotype matrix w.r.to
                    // the new ancient samples
                    if (!genetic_value_matrix.empty())
                        {
                            auto offset = individual * total_dim;
                            ancient_sample_genetic_value_matrix.insert(
                                end(ancient_sample_genetic_value_matrix),
                                begin(genetic_value_matrix) + offset,
                                begin(genetic_value_matrix) + offset + total_dim);
                        }
                }
        }

        /* This is the back end for importing mutations from tskit.
         * They have been decoded from metadata already.
         *
         * The procedure is:
         * * Add them to our table structures.
         * * Sort the tables.
         * * Traverse trees so that we can add mutations to genomes.
         *
         * Because the data come from a tskit TreeSequence, we
         * can make some simplifying assumptions:
         *
         * 1. Edges, etc., are already sorted.
         * 2. So we just need to sort mutations.
         */
        void
        set_mutations(const std::vector<Mutation> &mutations,
                      const std::vector<std::int32_t> &mutation_nodes)
        {
            if (this->is_simulating)
                {
                    throw std::runtime_error(
                        "cannot set mutations for a simulating population");
                }
            if (!this->mutations.empty())
                {
                    throw std::runtime_error("population has existing mutations");
                }
            this->mut_lookup.clear();
            this->tables->mutations.clear();
            this->tables->sites.clear();

            // Here, we cheat a bit
            this->haploid_genomes.clear();
            for (auto &dip : this->diploids)
                {
                    this->haploid_genomes.emplace_back(1);
                    dip.first = this->haploid_genomes.size() - 1;
                    this->haploid_genomes.emplace_back(1);
                    dip.second = this->haploid_genomes.size() - 1;
                }

            // Paranoia
            for (auto &g : this->haploid_genomes)
                {
                    if (g.n != 1)
                        {
                            throw std::runtime_error(
                                "all haploid_genomes must have a count of 1");
                        }
                }

            for (std::size_t i = 0; i < mutations.size(); ++i)
                {
                    if (this->mut_lookup.find(mutations[i].pos) != end(this->mut_lookup))
                        {
                            throw std::invalid_argument("duplicate mutation positions");
                        }
                    this->mutations.emplace_back(
                        Mutation(mutations[i].neutral, mutations[i].pos, mutations[i].s,
                                 mutations[i].h, -mutations[i].g, mutations[i].esizes,
                                 mutations[i].heffects, mutations[i].xtra));
                    this->mut_lookup.emplace(mutations[i].pos, i);
                    this->tables->emplace_back_site(mutations[i].pos, std::int8_t{0});
                    this->tables->emplace_back_mutation(
                        mutation_nodes[i], this->mutations.size() - 1,
                        this->tables->sites.size() - 1, std::int8_t{1},
                        mutations[i].neutral);
                }
            this->mcounts.resize(this->mutations.size());
            std::fill(begin(this->mcounts), end(this->mcounts), 0);
            fwdpp::ts::sort_mutation_table_and_rebuild_site_table(*this->tables);
            std::vector<fwdpp::ts::table_index_t> samples;
            std::unordered_map<fwdpp::ts::table_index_t, std::size_t> node_to_genome;
            for (std::size_t i = 0; i < this->diploid_metadata.size(); ++i)
                {
                    for (auto n : this->diploid_metadata[i].nodes)
                        {
                            if (node_to_genome.find(n) != end(node_to_genome))
                                {
                                    throw std::runtime_error(
                                        "node id use multiple times");
                                }
                            samples.push_back(n);
                        }
                    node_to_genome.emplace(this->diploid_metadata[i].nodes[0],
                                           this->diploids[i].first);
                    node_to_genome.emplace(this->diploid_metadata[i].nodes[1],
                                           diploids[i].second);
                }
            if (samples.empty())
                {
                    throw std::runtime_error(
                        "samples list for adding tskit mutations to genomes is empty");
                }

            auto sv = fwdpp::ts::site_visitor<fwdpp::ts::std_table_collection>(
                *this->tables, samples);
            auto site = sv();
            std::size_t nmutations_processed = 0;
            while (site != end(sv))
                {
                    auto mutations = sv.get_mutations();
                    if (mutations.second - mutations.first != 1)
                        {
                            throw std::runtime_error(
                                "expected exactly one mutation per site");
                        }
                    for (auto m = mutations.first; m != mutations.second; ++m)
                        {
                            fwdpp::ts::samples_iterator si(
                                sv.current_tree(), m->node,
                                fwdpp::ts::convert_sample_index_to_nodes(true));
                            auto s = si();
                            while (s != fwdpp::ts::NULL_INDEX)
                                {
                                    auto lookup = node_to_genome.find(s);
                                    if (lookup == end(node_to_genome))
                                        {
                                            throw std::runtime_error(
                                                "mutation could not be mapped to a "
                                                "sample node");
                                        }
                                    if (m->neutral)
                                        {
                                            this->haploid_genomes[lookup->second]
                                                .mutations.emplace_back(m->key);
                                        }
                                    else
                                        {
                                            this->haploid_genomes[lookup->second]
                                                .smutations.emplace_back(m->key);
                                        }
                                    this->mcounts[m->key] += 1;
                                    s = si();
                                }
                            nmutations_processed += 1;
                        }
                    site = sv();
                }
            if (nmutations_processed != this->mutations.size()
                || nmutations_processed != this->tables->mutations.size())
                {
                    throw std::runtime_error(
                        "failed to process the expected number of mutations");
                }
        }
    };
} // namespace fwdpy11
#endif
