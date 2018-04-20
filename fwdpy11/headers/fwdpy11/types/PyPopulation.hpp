#ifndef FWDPY11_PYPOPULATION_HPP__
#define FWDPY11_PYPOPULATION_HPP__

#include <pybind11/pybind11.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>

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

        virtual ~PyPopulation() = default;

        PyPopulation(PyPopulation &&) = default;
        PyPopulation(const PyPopulation &) = default;

        PyPopulation(fwdpp::uint_t N_)
            : fwdpp_base{ N_ }, N{ N_ }, generation{ 0 },
              popdata{ pybind11::none() }, popdata_user{ pybind11::none() }
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
              popdata_user{ pybind11::none() }
        {
        }

        virtual std::vector<std::size_t>
        add_mutations(typename fwdpp_base::mcont_t &mutations,
                      const std::vector<std::size_t> &individuals,
                      const std::vector<short> &gametes)
            = 0;
    };
}

#endif
