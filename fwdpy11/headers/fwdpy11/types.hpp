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
#include <fwdpy11/opaque/opaque_types.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpp/sugar.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <gsl/gsl_statistics_double.h>
#include <map>
#include <memory>
#include <vector>
#include <stdexcept>
#include <fwdpy11/serialization.hpp>

namespace fwdpy11
{
    //! Allows serialization of diploids.
    struct diploid_writer
    {
        using result_type = void;
        template <typename diploid_t, typename streamtype>
        inline result_type
        operator()(const diploid_t &dip, streamtype &o) const
        {
            KTfwd::fwdpp_internal::scalar_writer()(o, &dip.g);
            KTfwd::fwdpp_internal::scalar_writer()(o, &dip.e);
            KTfwd::fwdpp_internal::scalar_writer()(o, &dip.w);
            KTfwd::fwdpp_internal::scalar_writer()(o, &dip.label);
        }
    };

    //! Allows de-serialization of diploids.
    struct diploid_reader
    {
        using result_type = void;
        template <typename diploid_t, typename streamtype>
        inline result_type
        operator()(diploid_t &dip, streamtype &i) const
        {
            KTfwd::fwdpp_internal::scalar_reader()(i, &dip.g);
            KTfwd::fwdpp_internal::scalar_reader()(i, &dip.e);
            KTfwd::fwdpp_internal::scalar_reader()(i, &dip.w);
            KTfwd::fwdpp_internal::scalar_reader()(i, &dip.label);
        }
    };

    struct singlepop_t : public KTfwd::singlepop<KTfwd::popgenmut, diploid_t>
    /*!
      \brief Single-deme object where mutations have single effect size and
      dominance.

      This is the C++ representation of a single-deme simulation where a
      KTfwd::popgenmut
      is the mutation type, and a custom diploid (fwdpy11::diploid_t) is the
      diploid type.

      This type inherits from fwdpp's type KTfwd::singlepop, and the main
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

      \note Internally, fwdppy uses extinct elements in containers for "object
      recycling."

      Further:

      1. For each diploid, first and second are indexes into gametes
      2. Within each gamete, mutations and smutations contain indexes to
      "neutral" and "selected" mutations, resepectively, in mutations
    */
    {
        //! Typedef for base type
        using base = KTfwd::singlepop<KTfwd::popgenmut, diploid_t>;
        //! The current generation.  Start counting from zero
        unsigned generation;
        //! Constructor takes number of diploids as argument
        singlepop_t(const unsigned &N) : base(N), generation(0)
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

        singlepop_t(const std::string &s) : base(0) { this->deserialize(s); }

        singlepop_t(singlepop_t &&) = default;
        singlepop_t(const singlepop_t &) = default;
        singlepop_t &operator=(const singlepop_t &) = default;
        singlepop_t &operator=(singlepop_t &&) = default;

        std::string
        serialize() const
        {
            return serialization::serialize_details(
                this, KTfwd::mutation_writer(), fwdpy11::diploid_writer());
        }

        void
        deserialize(const std::string &s)
        {
            *this = serialization::deserialize_details<singlepop_t>()(
                s, KTfwd::mutation_reader<singlepop_t::mutation_t>(),
                fwdpy11::diploid_reader(), 1);
        }

        // int
        // tofile(const char *filename, bool append = false) const
        //{
        //    return fwdpy11::serialization::gzserialize_details(
        //        *this, KTfwd::mutation_writer(), fwdpy11::diploid_writer(),
        //        filename, append);
        //}

        // void
        // fromfile(const char *filename, std::size_t offset)
        //{
        //    *this = serialization::gzdeserialize_details<singlepop_t>()(
        //        KTfwd::mutation_reader<singlepop_t::mutation_t>(),
        //        fwdpy11::diploid_reader(), filename, offset, 0u);
        //}
    };

    struct metapop_t : public KTfwd::metapop<KTfwd::popgenmut, diploid_t>
    /*!
      \brief Multi-deme object where mutations have single effect size and
      dominance.

      This is the C++ representation of a single-deme simulation where a
      KTfwd::popgenmut
      is the mutation type, and a custom diploid (fwdpy11::diploid_t) is the
      diploid type.

      This type inherits from fwdpp's type KTfwd::singlepop, and the main
      documentation
      is found in the fwdpp reference manual.

      Here is a brief overview of the most important members of the base class.

      The format is name: type, comment

      mutations: vector<popgenmut>, contains a mix of extinct, segregating, and
      possibly fixed variants, depending on the simulation type.
      mcounts: vector<unsigned>, is the same length as mutations, and records
      the count (frequency, as an integer) of each mutation.
      gametes: vector<gamete>, contains a mix of extinct and extant gametes
      diploids: vector<vector<diploid>>, each vector is a deme
      fixations: vector<popgenmut>, is filled by simulations that "prune"
      fixations when they occur
      fixation_times: vector<unsigned>, is filled by simulations that "prune"
      fixations when they occur

      \note Internally, fwdppy uses extinct elements in containers for "object
      recycling."

      Further:

      1. For each diploid, first and second are indexes into gametes
      2. Within each gamete, mutations and smutations contain indexes to
      "neutral" and "selected" mutations, resepectively, in mutations
    */
    {
        //! Typedef for base type
        using base = KTfwd::metapop<KTfwd::popgenmut, diploid_t>;
        //! Current generation.  Start counting from 0
        unsigned generation;
        //! Constructor takes list of deme sizes are aregument
        explicit metapop_t(const std::vector<unsigned> &Ns)
            : base(&Ns[0], Ns.size()), generation(0)
        {
            // need to determine policy for how to label diploids :)
        }

        //! Constructor takes list of deme sizes are aregument
        explicit metapop_t(const std::initializer_list<unsigned> &Ns)
            : base(Ns), generation(0)
        {
        }

        //! Construct from a fwdpy11::singlepop_t
        explicit metapop_t(const singlepop_t &p)
            : base(p), generation(p.generation)
        {
        }
        // std::string
        // serialize() const
        //{
        //    return serialization::serialize_details(
        //        this, KTfwd::mutation_writer(), fwdpy11::diploid_writer());
        //}

        // void
        // deserialize(const std::string &s)
        //{
        //    *this = serialization::deserialize_details<metapop_t>()(
        //        s, KTfwd::mutation_reader<metapop_t::mutation_t>(),
        //        fwdpy11::diploid_reader(), std::vector<unsigned>(0u));
        //}

        // int
        // tofile(const char *filename, bool append = false) const
        //{
        //    return fwdpy11::serialization::gzserialize_details(
        //        *this, KTfwd::mutation_writer(), fwdpy11::diploid_writer(),
        //        filename, append);
        //}

        // void
        // fromfile(const char *filename, std::size_t offset)
        //{
        //    *this = serialization::gzdeserialize_details<metapop_t>()(
        //        KTfwd::mutation_reader<metapop_t::mutation_t>(),
        //        fwdpy11::diploid_reader(), filename, offset,
        //        std::vector<unsigned>(0u));
        //}
    };

    // Types based on KTfwd::generalmut_vec  //! Typedef for gamete type

    //! Typedef for gamete container
    using gcont_gm_vec_t = std::vector<gamete_t>;

    struct singlepop_gm_vec_t
        : public KTfwd::singlepop<KTfwd::generalmut_vec, diploid_t>
    /*!
      \brief Single-deme object where mutations contain vector<double> for
    internal data.
    ,
      See fwdpy11::singlepop_t documentation for details, which are the same as
    for this type.
    */
    {
        using base = KTfwd::singlepop<KTfwd::generalmut_vec, diploid_t>;
        unsigned generation;
        explicit singlepop_gm_vec_t(const unsigned &N) : base(N), generation(0)
        {
            if (!N)
                {
                    throw std::invalid_argument("population size must be > 0");
                }
        }
        explicit singlepop_gm_vec_t(const std::string &s) : base(0)
        {
            this->deserialize(s);
        }

        std::string
        serialize() const
        {
            return serialization::serialize_details(
                this, KTfwd::mutation_writer(), fwdpy11::diploid_writer());
        }

        void
        deserialize(const std::string &s)
        {
            *this = serialization::deserialize_details<singlepop_gm_vec_t>()(
                s, KTfwd::mutation_reader<singlepop_gm_vec_t::mutation_t>(),
                fwdpy11::diploid_reader(), 1u);
        }

        // int
        // tofile(const char *filename, bool append = false) const
        //{
        //    return fwdpy11::serialization::gzserialize_details(
        //        *this, KTfwd::mutation_writer(), fwdpy11::diploid_writer(),
        //        filename, append);
        //}

        // void
        // fromfile(const char *filename, std::size_t offset)
        //{
        //    *this = serialization::
        //        gzdeserialize_details<singlepop_gm_vec_t>()(
        //            KTfwd::mutation_reader<singlepop_gm_vec_t::mutation_t>(),
        //            fwdpy11::diploid_reader(), filename, offset, 0u);
        //}
    };

    // Types for multi-"locus" (multi-region) simulations
    using multilocus_diploid_t = std::vector<diploid_t>;

    // Have to use fwdpy11::diploid_t below, as GCC seems to get confused
    // otherwise...
    struct multilocus_t
        : public KTfwd::multiloc<KTfwd::popgenmut, fwdpy11::diploid_t>
    {
        using base = KTfwd::multiloc<KTfwd::popgenmut, fwdpy11::diploid_t>;
        unsigned generation, nloci;
        explicit multilocus_t(const unsigned N, const unsigned nloci_)
            : base(N, nloci_), generation(0), nloci(nloci_)
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
            : base(N, nloci_, locus_boundaries), generation(0), nloci(nloci_)
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

        explicit multilocus_t(const std::string &s) : base({ 0, 0 })
        {
            this->deserialize(s);
        }

        std::string
        serialize() const
        {
            return serialization::serialize_details(
                this, KTfwd::mutation_writer(), fwdpy11::diploid_writer());
        }

        void
        deserialize(const std::string &s)
        {
            *this = serialization::deserialize_details<multilocus_t>()(
                s, KTfwd::mutation_reader<multilocus_t::mutation_t>(),
                fwdpy11::diploid_reader(), 1, 1);
			if(!this->diploids.empty())
			{
				this->nloci=this->diploids[0].size();
			}
        }

        // int
        // tofile(const char *filename, bool append = false) const
        //{
        //    return fwdpy11::serialization::gzserialize_details(
        //        *this, KTfwd::mutation_writer(), fwdpy11::diploid_writer(),
        //        filename, append);
        //}

        // void
        // fromfile(const char *filename, std::size_t offset)
        //{
        //    *this = serialization::gzdeserialize_details<multilocus_t>()(
        //        KTfwd::mutation_reader<multilocus_t::mutation_t>(),
        //        fwdpy11::diploid_reader(), filename, offset, 0u, 0u);
        //}
    };
}

#endif
