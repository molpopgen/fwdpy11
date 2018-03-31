//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
/*!
  \file types.hpp

  \brief Wrappers around fwdpp types for the fwdpy Python package.
*/

#ifndef FWDPY11_TYPES__
#define FWDPY11_TYPES__

#include <pybind11/pybind11.h>
//#include <fwdpy11/opaque/opaque_types.hpp>
#include "types/typedefs.hpp"
#include "types/diploid.hpp"
#include <fwdpy11/rng.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
//#include <fwdpp/sugar/generalmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/mlocuspop.hpp>
#include <vector>
#include <stdexcept>
#include <fwdpy11/serialization.hpp>

namespace fwdpy11
{
    template <typename poptype, typename diploids_input,
              typename gametes_input, typename mutations_input>
    inline poptype
    create_wrapper(diploids_input &&diploids, gametes_input &&gametes,
                   mutations_input &&mutations)
    {
        return poptype(std::forward<diploids_input>(diploids),
                       std::forward<gametes_input>(gametes),
                       std::forward<mutations_input>(mutations));
    }

    template <typename poptype, typename diploids_input,
              typename gametes_input, typename mutations_input>
    inline poptype
    create_wrapper(diploids_input &&diploids, gametes_input &&gametes,
                   mutations_input &&mutations, mutations_input &fixations,
                   std::vector<fwdpp::uint_t> &fixation_times,
                   fwdpp::uint_t generation)
    {
        if (fixation_times.size() != fixations.size())
            {
                throw pybind11::value_error(
                    "length of fixation_times != length of fixations");
            }
        auto rv = create_wrapper<poptype>(
            std::forward<diploids_input>(diploids),
            std::forward<gametes_input>(gametes),
            std::forward<mutations_input>(mutations));
        rv.fixations.swap(fixations);
        rv.fixation_times.swap(fixation_times);
        rv.generation = generation;
        return rv;
    }

    //! Allows serialization of diploids.
    //template <int VERSION> struct diploid_writer
    //{
    //    using result_type = void;
    //    // This should really be constexpr.
    //    // Figure it out later:
    //    const int v;
    //    diploid_writer() : v(VERSION) {}
    //    template <typename diploid_t, typename streamtype>
    //    inline result_type
    //    operator()(const diploid_t &dip, streamtype &o) const
    //    {
    //        fwdpp::fwdpp_internal::scalar_writer()(o, &dip.g);
    //        fwdpp::fwdpp_internal::scalar_writer()(o, &dip.e);
    //        fwdpp::fwdpp_internal::scalar_writer()(o, &dip.w);
    //        fwdpp::fwdpp_internal::scalar_writer()(o, &dip.label);
    //    }
    //};

    //! Allows de-serialization of diploids.
    // struct diploid_reader
    // {
    //     using result_type = void;
    //     template <typename diploid_t, typename streamtype>
    //     inline result_type
    //     operator()(diploid_t &dip, streamtype &i, const int version) const
    //     {
    //         fwdpp::fwdpp_internal::scalar_reader()(i, &dip.g);
    //         fwdpp::fwdpp_internal::scalar_reader()(i, &dip.e);
    //         fwdpp::fwdpp_internal::scalar_reader()(i, &dip.w);
    //         fwdpp::fwdpp_internal::scalar_reader()(i, &dip.label);
    //     }
    // };

    struct singlepop_t : public fwdpp::slocuspop<fwdpp::popgenmut, Diploid>
    /*!
      \brief Single-deme object where mutations have single effect size and
      dominance.

      This is the C++ representation of a single-deme simulation where a
      fwdpp::popgenmut
      is the mutation type, and a custom diploid (fwdpy11::Diploid) is the
      diploid type.

      This type inherits from fwdpp's type fwdpp::singlepop, and the main
      documentation
      is found in the fwdpp reference manual.

      Here is a brief overview of the most important members of the base class.

      The format is name: type, comment

      mutations: vector<popgenmut>, contains a mix of extinct, segregating, and
      possibly fixed variants, depending on the simulation type.
      mcounts: vector<unsigned>, is the same length as mutations, and records
      the count (frequency, as an integer) of each mutation.
      gametes: vector<gamete>, contains a mix of extinct and extant gametes
      diploids: vector<diploid>
      fixations: vector<popgenmut>, is filled by simulations that "prune"
      fixations when they occur
      fixation_times: vector<unsigned>, is filled by simulations that "prune"
      fixations when they occur

      \note Internally, fwdpp uses extinct elements in containers for "object
      recycling."

      Further:

      1. For each diploid, first and second are indexes into gametes
      2. Within each gamete, mutations and smutations contain indexes to
      "neutral" and "selected" mutations, resepectively, in mutations
    */
    {
        //! Typedef for base type
        using base = fwdpp::slocuspop<fwdpp::popgenmut, Diploid>;
        //! The current generation.  Start counting from zero
        unsigned generation;
        //! A Python object that sims may (optionally) use
        pybind11::object popdata;
        //! A Python objeft that users may access during a simulation
        pybind11::object popdata_user;
        //! Constructor takes number of diploids as argument
        singlepop_t(const unsigned &N)
            : base(N),
              generation(0), popdata{ pybind11::none() }, popdata_user{
                  pybind11::none()
              }
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

        // Perfect-forwarding constructor:
        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        singlepop_t(diploids_input &&diploids, gametes_input &&gametes,
                    mutations_input &&mutations)
            : base(std::forward<diploids_input>(diploids),
                   std::forward<gametes_input>(gametes),
                   std::forward<mutations_input>(mutations)),
              generation{ 0 }, popdata{ pybind11::none() }, popdata_user{
                  pybind11::none()
              }
        {
        }

        singlepop_t(const std::string &s) : base(0) { this->deserialize(s); }

        singlepop_t(singlepop_t &&) = default;
        singlepop_t(const singlepop_t &) = default;
        singlepop_t &operator=(const singlepop_t &) = default;
        singlepop_t &operator=(singlepop_t &&) = default;

        static singlepop_t
        create(base::dipvector_t &diploids, base::gcont_t &gametes,
               base::mcont_t &mutations)
        {
            return create_wrapper<singlepop_t>(
                std::move(diploids), std::move(gametes), std::move(mutations));
        }

        static singlepop_t
        create_with_fixations(base::dipvector_t &diploids,
                              base::gcont_t &gametes, base::mcont_t &mutations,
                              base::mcont_t &fixations,
                              std::vector<fwdpp::uint_t> &fixation_times,
                              const fwdpp::uint_t generation)
        {
            return create_wrapper<singlepop_t>(
                std::move(diploids), std::move(gametes), std::move(mutations),
                fixations, fixation_times, generation);
        }

        std::string
        serialize() const
        {
            return serialization::serialize_details(this);
        }

        void
        deserialize(const std::string &s)
        {
            *this = serialization::deserialize_details<singlepop_t>()(s,1);
        }
    };

    // Types based on fwdpp::generalmut_vec  //! Typedef for gamete type

    //! Typedef for gamete container
    // using gcont_gm_vec_t = std::vector<gamete_t>;

    // struct singlepop_gm_vec_t
    //     : public fwdpp::singlepop<fwdpp::generalmut_vec, diploid_t>
    // /*!
    //   \brief Single-deme object where mutations contain vector<double> for
    // internal data.
    // ,
    //   See fwdpy11::singlepop_t documentation for details, which are the
    // same as
    // for this type.
    // */
    // {
    //     using base = fwdpp::singlepop<fwdpp::generalmut_vec, diploid_t>;
    //     unsigned generation;
    //     pybind11::object popdata;
    //     //! A Python objeft that users may access during a simulation
    //     pybind11::object popdata_user;
    //     //! Constructor takes number of diploids as argument
    //     explicit singlepop_gm_vec_t(const unsigned &N)
    //         : base(N), generation(0), popdata{ pybind11::none() },
    //           popdata_user{ pybind11::none() }
    //     {
    //         if (!N)
    //             {
    //                 throw std::invalid_argument("population size must be >
    //                 0");
    //             }
    //     }

    //     // Perfect-forwarding constructor:
    //     template <typename diploids_input, typename gametes_input,
    //               typename mutations_input>
    //     explicit singlepop_gm_vec_t(diploids_input &&diploids,
    //                                 gametes_input &&gametes,
    //                                 mutations_input &&mutations)
    //         : base(std::forward<diploids_input>(diploids),
    //                std::forward<gametes_input>(gametes),
    //                std::forward<mutations_input>(mutations)),
    //           generation{ 0 }, popdata{ pybind11::none() },
    //           popdata_user{ pybind11::none() }
    //     {
    //     }

    //     explicit singlepop_gm_vec_t(const std::string &s)
    //         : base(0), popdata{ pybind11::none() },
    //           popdata_user{ pybind11::none() }
    //     {
    //         this->deserialize(s);
    //     }

    //     static singlepop_gm_vec_t
    //     create(base::dipvector_t &diploids, base::gcont_t &gametes,
    //            base::mcont_t &mutations)
    //     {
    //         return create_wrapper<singlepop_gm_vec_t>(
    //             std::move(diploids), std::move(gametes),
    //             std::move(mutations));
    //     }

    //     static singlepop_gm_vec_t
    //     create_with_fixations(base::dipvector_t &diploids,
    //                           base::gcont_t &gametes, base::mcont_t
    //                           &mutations,
    //                           base::mcont_t &fixations,
    //                           std::vector<fwdpp::uint_t> &fixation_times,
    //                           const fwdpp::uint_t generation)
    //     {
    //         return create_wrapper<singlepop_gm_vec_t>(
    //             std::move(diploids), std::move(gametes),
    //             std::move(mutations),
    //             fixations, fixation_times, generation);
    //     }
    //     std::string
    //     serialize() const
    //     {
    //         return serialization::serialize_details(
    //             this, fwdpp::mutation_writer(),
    //             fwdpy11::diploid_writer<fwdpy11::serialization::magic()>());
    //     }

    //     void
    //     deserialize(const std::string &s)
    //     {
    //         *this =
    //         serialization::deserialize_details<singlepop_gm_vec_t>()(
    //             s, fwdpp::mutation_reader<singlepop_gm_vec_t::mutation_t>(),
    //             fwdpy11::diploid_reader(), 1u);
    //     }

    //     // int
    //     // tofile(const char *filename, bool append = false) const
    //     //{
    //     //    return fwdpy11::serialization::gzserialize_details(
    //     //        *this, fwdpp::mutation_writer(),
    //     //        fwdpy11::diploid_writer(),
    //     //        filename, append);
    //     //}

    //     // void
    //     // fromfile(const char *filename, std::size_t offset)
    //     //{
    //     //    *this = serialization::
    //     //        gzdeserialize_details<singlepop_gm_vec_t>()(
    //     // fwdpp::mutation_reader<singlepop_gm_vec_t::mutation_t>(),
    //     //            fwdpy11::diploid_reader(), filename, offset, 0u);
    //     //}
    // };

    // Types for multi-"locus" (multi-region) simulations
    using multilocus_Diploid = std::vector<Diploid>;

    // Have to use fwdpy11::diploid_t below, as GCC seems to get confused
    // otherwise...
    struct multilocus_t
        : public fwdpp::mlocuspop<fwdpp::popgenmut, fwdpy11::Diploid>
    {
        using base = fwdpp::mlocuspop<fwdpp::popgenmut, fwdpy11::Diploid>;
        unsigned generation, nloci;
        pybind11::object popdata;
        //! A Python objeft that users may access during a simulation
        pybind11::object popdata_user;
        //! Constructor takes number of diploids as argument
        explicit multilocus_t(const unsigned N, const unsigned nloci_)
            : base(N, nloci_), generation(0),
              nloci(nloci_), popdata{ pybind11::none() }, popdata_user{
                  pybind11::none()
              }
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
        explicit multilocus_t(
            const unsigned N, const unsigned nloci_,
            const std::vector<std::pair<double, double>> &locus_boundaries)
            : base(N, nloci_, locus_boundaries), generation(0),
              nloci(nloci_), popdata{ pybind11::none() }, popdata_user{
                  pybind11::none()
              }
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

        // Perfect-forwarding constructor:
        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        multilocus_t(diploids_input &&diploids, gametes_input &&gametes,
                     mutations_input &&mutations)
            : base(std::forward<diploids_input>(diploids),
                   std::forward<gametes_input>(gametes),
                   std::forward<mutations_input>(mutations)),
              generation{ 0 }, popdata{ pybind11::none() }, popdata_user{
                  pybind11::none()
              }
        {
        }

        explicit multilocus_t(const std::string &s) : base({ 0, 0 })
        {
            this->deserialize(s);
        }

        static multilocus_t
        create(base::dipvector_t &diploids, base::gcont_t &gametes,
               base::mcont_t &mutations)
        {
            return create_wrapper<multilocus_t>(
                std::move(diploids), std::move(gametes), std::move(mutations));
        }

        static multilocus_t
        create_with_fixations(base::dipvector_t &diploids,
                              base::gcont_t &gametes, base::mcont_t &mutations,
                              base::mcont_t &fixations,
                              std::vector<fwdpp::uint_t> &fixation_times,
                              const fwdpp::uint_t generation)
        {
            return create_wrapper<multilocus_t>(
                std::move(diploids), std::move(gametes), std::move(mutations),
                fixations, fixation_times, generation);
        }
        std::string
        serialize() const
        {
            return serialization::serialize_details(this);
        }

        void
        deserialize(const std::string &s)
        {
            *this
                = serialization::deserialize_details<multilocus_t>()(s, 1, 1);
            if (!this->diploids.empty())
                {
                    this->nloci = this->diploids[0].size();
                }
        }

        // int
        // tofile(const char *filename, bool append = false) const
        //{
        //    return fwdpy11::serialization::gzserialize_details(
        //        *this, fwdpp::mutation_writer(),
        //        fwdpy11::diploid_writer(),
        //        filename, append);
        //}

        // void
        // fromfile(const char *filename, std::size_t offset)
        //{
        //    *this = serialization::gzdeserialize_details<multilocus_t>()(
        //        fwdpp::mutation_reader<multilocus_t::mutation_t>(),
        //        fwdpy11::diploid_reader(), filename, offset, 0u, 0u);
        //}
    };
}

#endif
