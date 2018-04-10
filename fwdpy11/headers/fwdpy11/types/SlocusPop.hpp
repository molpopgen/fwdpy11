#ifndef FWDPY11_SLOCUSPOP_HPP__
#define FWDPY11_SLOCUSPOP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
#include "create_pops.hpp"
#include <stdexcept>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace fwdpy11
{
    class SlocusPop : public Population
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
        using dipvector_t = std::vector<Diploid>;
        using diploid_t = dipvector_t::value_type;
        using popbase_t = Population;
        using popmodel_t = fwdpp::sugar::SINGLELOC_TAG;
        using fitness_t
            = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                           typename popbase_t::gcont_t,
                                           typename popbase_t::mcont_t>;

        dipvector_t diploids;

        SlocusPop(SlocusPop &&) = default;
        SlocusPop(const SlocusPop &) = default;
        SlocusPop &operator=(SlocusPop &&) = default;
        SlocusPop &operator=(const SlocusPop &) = default;

        // Constructors for Python
        SlocusPop(const fwdpp::uint_t N)
            : Population{ N }, diploids(N, { 0, 0 })
        {
            if (!N)
                {
                    throw std::invalid_argument("population size must be > 0");
                }
            std::size_t label = 0;
            for (auto &&d : this->diploids)
                {
                    d.label = label++;
                }
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        explicit SlocusPop(diploids_input &&d, gametes_input &&g,
                           mutations_input &&m)
            : Population(static_cast<fwdpp::uint_t>(d.size()),
                         std::forward<gametes_input>(g),
                         std::forward<mutations_input>(m), 100),
              diploids(std::forward<diploids_input>(d))
        //! Constructor for pre-determined population status
        {
            this->process_individual_input();
        }

        bool
        operator==(const SlocusPop &rhs) const
        {
            return this->diploids == rhs.diploids && popbase_t::is_equal(rhs);
        };

        void
        clear()
        {
            diploids.clear();
            popbase_t::clear_containers();
        }

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &gametes)
        {
            return {};
        }
    };
}
#endif
