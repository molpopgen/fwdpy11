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
#ifndef FWDPY11_SERIALIZATION_DIPLOID_HPP__
#define FWDPY11_SERIALIZATION_DIPLOID_HPP__

#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/io/diploid.hpp>
#include <fwdpp/io/scalar_serialization.hpp>

// Make Diploid compatible with fwdpp's
// serialization API:

namespace fwdpp
{
    namespace io
    {
         template <> struct serialize_diploid<fwdpy11::DiploidGenotype>
         {
             template <typename streamtype>
             inline void
             operator()(streamtype& buffer, const fwdpy11::DiploidGenotype& dip) const
             {
                 fwdpp::io::scalar_writer w;
                 w(buffer, &dip.first);
                 w(buffer, &dip.second);
             }
         };

         template <> struct deserialize_diploid<fwdpy11::DiploidGenotype>
         {
             template <typename streamtype>
             inline void
             operator()(streamtype& buffer, fwdpy11::DiploidGenotype& dip) const
             {
                 fwdpp::io::scalar_reader r;
                 r(buffer, &dip.first);
                 r(buffer, &dip.second);
             }
         };
    }
}

#endif

