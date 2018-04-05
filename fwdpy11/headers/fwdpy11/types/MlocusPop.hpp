#ifndef FWDPY11_MLOCUSPOP_HPP__
#define FWDPY11_MLOCUSPOP_HPP__

#include "Population.hpp"
#include "Diploid.hpp"
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
        MlocusPop(const fwdpp::uint_t N, fwdpp::uint_t nloci_)
            : Population{ N },
              diploids(N, diploid_t(nloci_, Diploid{ 0, 0 })), nloci{ nloci_ },
              locus_boundaries{}
        {
        }

        MlocusPop(
            const fwdpp::uint_t N, fwdpp::uint_t nloci_,
            const std::vector<std::pair<double, double>> &locus_boundaries_)
            : Population{ N },
              diploids(N, diploid_t(nloci_, Diploid{ 0, 0 })), nloci{ nloci_ },
              locus_boundaries{ locus_boundaries_ }
        {
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        explicit MlocusPop(diploids_input &&d, gametes_input &&g,
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
        operator==(const MlocusPop &rhs) const
        {
            return this->diploids == rhs.diploids && popbase_t::is_equal(rhs);
        };

        void
        clear()
        {
            diploids.clear();
            popbase_t::clear_containers();
        }
    };
}
#endif
