#ifndef FWDPY11_MLOCUSPOP_HPP__
#define FWDPY11_MLOCUSPOP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <fwdpp/sugar/add_mutation.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace fwdpy11
{
    class MlocusPop : public Population
    {
      private:
        void
        process_individual_input()
        {
            std::vector<fwdpp::uint_t> gcounts(this->gametes.size(), 0);
            for (auto &&dip : diploids)
                {
                    for (auto &&locus : dip)
                        {
                            this->validate_individual_keys(locus.first);
                            this->validate_individual_keys(locus.second);
                            gcounts[locus.first]++;
                            gcounts[locus.second]++;
                        }
                }
            this->validate_gamete_counts(gcounts);
        }

      public:
        using dipvector_t = std::vector<std::vector<Diploid>>;
        using diploid_t = dipvector_t::value_type;
        using popbase_t = Population;
        using popmodel_t = fwdpp::sugar::MULTILOC_TAG;
        using fitness_t
            = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                           typename popbase_t::gcont_t,
                                           typename popbase_t::mcont_t>;

        dipvector_t diploids;
        fwdpp::uint_t nloci;
        std::vector<std::pair<double, double>> locus_boundaries;

        MlocusPop(MlocusPop &&) = default;
        MlocusPop(const MlocusPop &) = default;
        MlocusPop &operator=(MlocusPop &&) = default;
        MlocusPop &operator=(const MlocusPop &) = default;

        // Constructors for Python
        MlocusPop(const fwdpp::uint_t N,
                  std::vector<std::pair<double, double>> locus_boundaries_)
            : Population{ N }, diploids(N, diploid_t(locus_boundaries_.size(),
                                                     Diploid{ 0, 0 })),
              nloci{ static_cast<fwdpp::uint_t>(locus_boundaries_.size()) },
              locus_boundaries{ std::move(locus_boundaries_) }
        {
            if (!N)
                {
                    throw std::invalid_argument("population size must be > 0");
                }
            if (!nloci)
                {
                    throw std::invalid_argument("number of loci must be > 0");
                }
            std::size_t label = 0;
            for (auto &&d : this->diploids)
                {
                    d[0].label = label++;
                }
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        explicit MlocusPop(
            diploids_input &&d, gametes_input &&g, mutations_input &&m,
            std::vector<std::pair<double, double>> locus_boundaries_)
            : Population(static_cast<fwdpp::uint_t>(d.size()),
                         std::forward<gametes_input>(g),
                         std::forward<mutations_input>(m), 100),
              diploids(std::forward<diploids_input>(d)),
              nloci{ static_cast<fwdpp::uint_t>(locus_boundaries_.size()) }, locus_boundaries{ std::move(
                                                     locus_boundaries_) }
        //! Constructor for pre-determined population status
        {
            if (diploids.at(0).size() != nloci)
                {
                    throw std::invalid_argument("diploid genotypes "
                                                "inconsistent wit number of "
                                                "locus boundaries");
                }
            this->process_individual_input();
        }

        bool
        operator==(const MlocusPop &rhs) const
        {
            return this->diploids == rhs.diploids
                   && this->locus_boundaries == rhs.locus_boundaries
                   && popbase_t::is_equal(rhs);
        };

        void
        clear()
        {
            diploids.clear();
            locus_boundaries.clear();
            popbase_t::clear_containers();
        }

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &new_mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &gametes)
        {
            std::unordered_set<double> poschecker;
            std::vector<std::size_t> loci;
            for (const auto &m : new_mutations)
                {
                    auto x = std::lower_bound(
                        std::begin(locus_boundaries),
                        std::end(locus_boundaries), m.pos,
                        [](const std::pair<double, double> &p,
                           const double value) { return p.second < value; });
                    if (x == std::end(locus_boundaries))
                        {
                            throw std::invalid_argument("mutation position "
                                                        "not within locus "
                                                        "boundaries");
                        }
                    if (m.pos >= x->second)
                        {
                            throw std::invalid_argument(
                                "mutation position greater than right-hand "
                                "locus boundary");
                        }
                    loci.push_back(x - std::begin(locus_boundaries));
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

            std::size_t locus = 0;
            for (auto &i : new_mutations)
                {
                    auto pos = i.pos;
                    // remaining preconditions get checked by fwdpp:
                    rv.push_back(fwdpp::add_mutation(*this, loci[locus++],
                                                     individuals, gametes,
                                                     std::move(i)));
                    this->mut_lookup.insert(pos);
                }
            return rv;
        }
    };
}
#endif
