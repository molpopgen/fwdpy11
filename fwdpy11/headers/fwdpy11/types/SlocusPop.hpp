#ifndef FWDPY11_SLOCUSPOP_HPP__
#define FWDPY11_SLOCUSPOP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
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

        // static SlocusPop
        // create(dipvector_t &diploids, gcont_t &gametes, mcont_t &mutations)
        //{
        //    return create_wrapper<SlocusPop>(
        //        std::move(diploids), std::move(gametes),
        //        std::move(mutations));
        //}

        // static SlocusPop
        // create_with_fixations(dipvector_t &diploids, gcont_t &gametes,
        //                      mcont_t &mutations, mcont_t &fixations,
        //                      std::vector<fwdpp::uint_t> &fixation_times,
        //                      const fwdpp::uint_t generation)
        //{
        //    return create_wrapper<SlocusPop>(
        //        std::move(diploids), std::move(gametes),
        //        std::move(mutations),
        //        fixations, fixation_times, generation);
        //}
    };
}
#endif
