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
#ifndef FWDPY11_MUTATION_TYPE_HPP__
#define FWDPY11_MUTATION_TYPE_HPP__

/*
 * This file started off via a copy
 * of fwdpp's popgenmut.hpp
 * by Kevin Thornton.
 *
 * This file defines the minimal C++ API
 * for the type.
 *
 * To serialize instaces of this type
 * using fwdpp's machinery:
 *
 * #include <fwdpy11/serialization/Mutation.hpp>
 */

#include <tuple>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>

namespace fwdpy11
{
    struct Mutation : public fwdpp::mutation_base
    ///! The fwdpy11 mutation type
    {
        //! The generation when the mutation arose
        fwdpp::uint_t g;
        //! Effect size.  We call it 's' so
        // that we can use fwdpp's genetic value toolkit
        double s;
        //! Dominance of the mutation
        double h;
        std::vector<double> esizes, heffects;
        //! Alias for tuple type that can be used for object construction
        using constructor_tuple
            = std::tuple<double, double, double, unsigned, std::uint16_t>;
        //! Alias for tuple type that accommodates variable effect sizes
        using constructor_tuple_variable_effects
            = std::tuple<double, double, double, unsigned, std::vector<double>,
                         std::vector<double>, std::uint16_t>;

        /*!
          Constructor for constant effect size sims.

          \param pos_ Mutation position
          \param s_ Selection coefficient
          \param h_ Dominance coefficient
          \param g_ Generation when mutation arose
          \param x_ Value to assign to mutation_base::xtra
        */
        Mutation(const double &pos_, const double &s_, const double &h_,
                 const unsigned &g_, const std::uint16_t x_ = 0) noexcept
            : mutation_base(pos_, (s_ == 0.) ? true : false, x_), g(g_), s(s_),
              h(h_), esizes{}, heffects{}
        {
        }

        /*!
          Constructor for constant effect size + variable effect size sims.

          \param pos_ Mutation position
          \param s_ Selection coefficient
          \param h_ Dominance coefficient
          \param g_ Generation when mutation arose
          \param x_ Value to assign to mutation_base::xtra
          \param esizes_ Vector of effect sizes
          \param heffects_ Vector of heterozygous effects
        */
        template <typename vectype>
        Mutation(const double &pos_, const double &s_, const double &h_,
                 const unsigned &g_, vectype &&esizes_, vectype &&heffects_,
                 const std::uint16_t x_ = 0) noexcept
            : fwdpp::mutation_base(pos_, true, x_), g(g_), s(s_), h(h_),
              esizes(std::forward<vectype>(esizes_)),
              heffects(std::forward<vectype>(heffects_))
        {
            this->neutral
                = ((s == 0.0)
                   && std::all_of(std::begin(this->esizes),
                                  std::end(this->esizes),
                                  [](const double d) { return d == 0.; }));
        }

        Mutation(constructor_tuple t) noexcept
            : mutation_base(std::get<0>(t),
                            (std::get<1>(t) == 0.) ? true : false,
                            std::get<4>(t)),
              g(std::get<3>(t)), s(std::get<1>(t)),
              h(std::get<2>(t)), esizes{}, heffects{}
        {
            this->neutral
                = ((s == 0.0)
                   && std::all_of(std::begin(this->esizes),
                                  std::end(this->esizes),
                                  [](const double d) { return d == 0.; }));
        }

        Mutation(constructor_tuple_variable_effects t) noexcept
            : mutation_base(std::get<0>(t),
                            (std::get<1>(t) == 0.) ? true : false,
                            std::get<6>(t)),
              g(std::get<3>(t)), s(std::get<1>(t)), h(std::get<2>(t)),
              esizes(std::get<4>(t)), heffects(std::get<5>(t))
        {
            this->neutral
                = ((s == 0.0)
                   && std::all_of(std::begin(this->esizes),
                                  std::end(this->esizes),
                                  [](const double d) { return d == 0.; }));
        }

        bool
        operator==(const Mutation &rhs) const
        {
            return std::tie(this->g, this->s, this->h, this->esizes,
                            this->heffects)
                       == std::tie(rhs.g, rhs.s, rhs.h, rhs.esizes,
                                   rhs.heffects)
                   && is_equal(rhs);
        }
    };
}

#endif
