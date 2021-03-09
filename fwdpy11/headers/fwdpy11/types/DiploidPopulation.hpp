#ifndef FWDPY11_POP_HPP__
#define FWDPY11_POP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <fwdpp/poptypes/tags.hpp>

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
            const auto update = [&unodes,&n](auto i)
            {
                if(unodes.find(i) == end(unodes))
                {
                    n.push_back(i);
                    unodes.insert(i);
                }
            };
            n.clear();
            for (auto &md : metadata)
                {
                    for(auto i : md.nodes)
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
            if (this->tables->genome_length()
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
                                    tables->nodes[md.nodes[0]].deme
                                        = deme_label;
                                    tables->nodes[md.nodes[1]].deme
                                        = deme_label;
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
            : Population{ N, length }, diploids(N, { 0, 0 }),
              diploid_metadata(N), ancient_sample_metadata{}
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
    };
} // namespace fwdpy11
#endif
