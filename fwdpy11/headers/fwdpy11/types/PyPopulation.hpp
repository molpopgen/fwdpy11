#ifndef FWDPY11_PYPOPULATION_HPP__
#define FWDPY11_PYPOPULATION_HPP__

#include <tuple>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
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

      public:
        using fwdpp_base
            = fwdpp::sugar::popbase<mutation_type, mcont, gcont, mvector,
                                    ftvector, lookup_table_type>;
        fwdpp::uint_t N;
        fwdpp::uint_t generation;
        pybind11::object popdata;
        pybind11::object popdata_user;

        std::vector<DiploidMetadata> diploid_metadata;

        virtual ~PyPopulation() = default;

        PyPopulation(PyPopulation &&) = default;
        PyPopulation(const PyPopulation &) = default;

        PyPopulation(fwdpp::uint_t N_)
            : fwdpp_base{ N_ }, N{ N_ }, generation{ 0 },
              popdata{ pybind11::none() }, popdata_user{ pybind11::none() },
              diploid_metadata(N)
        {
        }

        template <typename gametes_input, typename mutations_input>
        explicit PyPopulation(
            const fwdpp::uint_t N_, gametes_input &&g, mutations_input &&m,
            typename fwdpp_base::gamete_t::mutation_container::size_type
                reserve_size)
            : fwdpp_base{ std::forward<gametes_input>(g),
                          std::forward<mutations_input>(m), reserve_size },
              N{ N_ }, generation{ 0 }, popdata{ pybind11::none() },
              popdata_user{ pybind11::none() },
              diploid_metadata(N)
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
    };
} // namespace fwdpy11

#endif
