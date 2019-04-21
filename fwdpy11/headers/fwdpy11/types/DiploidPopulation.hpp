#ifndef FWDPY11_POP_HPP__
#define FWDPY11_POP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include "create_pops.hpp"
#include <stdexcept>
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
            std::vector<fwdpp::uint_t> gcounts(this->gametes.size(), 0);
            for (auto &&dip : diploids)
                {
                    this->validate_individual_keys(dip.first);
                    this->validate_individual_keys(dip.second);
                    gcounts[dip.first]++;
                    gcounts[dip.second]++;
                }
            this->validate_gamete_counts(gcounts);
        }

      public:
        using dipvector_t = std::vector<DiploidGenotype>;
        using diploid_t = dipvector_t::value_type;
        using popbase_t = Population;
        using popmodel_t = fwdpp::poptypes::SINGLELOC_TAG;
        using fitness_t
            = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                           typename popbase_t::gcont_t,
                                           typename popbase_t::mcont_t>;

        dipvector_t diploids;

        // Constructors for Python
        DiploidPopulation(const fwdpp::uint_t N, const double length)
            : Population{ N, length }, diploids(N, { 0, 0 })
        {
            if (!N)
                {
                    throw std::invalid_argument("population size must be > 0");
                }
            std::size_t label = 0;
            for (auto &&d : this->diploid_metadata)
                {
                    d.label = label++;
                    d.w = 1.0;
                    if (length == std::numeric_limits<double>::max())
                        {
                            d.nodes[0] = d.nodes[1] = -1;
                        }
                }
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        explicit DiploidPopulation(diploids_input &&d, gametes_input &&g,
                           mutations_input &&m)
            : Population(static_cast<fwdpp::uint_t>(d.size()),
                         std::forward<gametes_input>(g),
                         std::forward<mutations_input>(m), 100),
              diploids(std::forward<diploids_input>(d))
        //! Constructor for pre-determined population status
        {
            this->process_individual_input();
        }

        ~DiploidPopulation() = default;
        DiploidPopulation(DiploidPopulation&&)=default;
        DiploidPopulation(const DiploidPopulation&)=default;
        DiploidPopulation&operator=(const DiploidPopulation&)=default;
        DiploidPopulation&operator=(DiploidPopulation&&)=default;

        bool
        operator==(const DiploidPopulation &rhs) const
        {
            return this->diploids == rhs.diploids && popbase_t::is_equal(rhs)
                   && popbase_t::tables_equal(rhs);
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
                      const std::vector<short> &gametes)
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
                    auto idx = fwdpp::add_mutation(*this, individuals, gametes,
                                                   std::move(i));

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
    };
} // namespace fwdpy11
#endif
