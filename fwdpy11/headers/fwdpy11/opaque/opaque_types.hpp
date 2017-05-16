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
#ifndef FWDPY11_OPAQUE_TYPES_HPP__
#define FWDPY11_OPAQUE_TYPES_HPP__

#include <pybind11/stl.h>
#include <cstdint>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>

namespace fwdpy11
{
    //! Typedef for mutation container
    using mcont_t = std::vector<KTfwd::popgenmut>;
    //! Typedef for gamete type
    using gamete_t = KTfwd::gamete;
    //! Typedef for gamete container
    using gcont_t = std::vector<gamete_t>;

    struct diploid_t
    /*!
      \brief Custom diploid type.
    */
    {
        using first_type = std::size_t;
        using second_type = std::size_t;
        //! First gamete.  A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        first_type first;
        //! Second gamete. A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        second_type second;
        //! 64 bits of data to do stuff with.  Initialized to zero upon
        //! construction
        std::size_t label;
        //! Genetic component of trait value.  This is not necessarily written
        //! to by a simulation.
        double g;
        //! Random component of trait value.  This is not necessarily written
        //! to by a simulation.
        double e;
        //! Fitness.  This is not necessarily written to by a simulation.
        double w;
        //! Constructor
        diploid_t() noexcept : first(first_type()),
                               second(second_type()),
                               label(0),
                               g(0.),
                               e(0.),
                               w(1.)
        {
        }
        //! Construct from two indexes to gametes
        diploid_t(first_type g1, first_type g2) noexcept : first(g1),
                                                           second(g2),
                                                           label(0),
                                                           g(0.),
                                                           e(0.),
                                                           w(1.)
        {
        }

        inline bool
        operator==(const diploid_t& dip) const noexcept
        //! Required for py::bind_vector
        {
            return this->first == dip.first && this->second == dip.second
                   && this->w == dip.w && this->g == dip.g && this->e == dip.e
                   && this->label == dip.label;
        }
    };

    //! Typedef for container of diploids
    using dipvector_t = std::vector<diploid_t>;
}

PYBIND11_MAKE_OPAQUE(fwdpy11::dipvector_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::mcont_t);
PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::generalmut_vec>);

#endif
