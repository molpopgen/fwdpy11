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
#ifndef FWDPY11_SERIALIZATION_MUTATION_HPP__
#define FWDPY11_SERIALIZATION_MUTATION_HPP__

#include <fwdpy11/types/Mutation.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/scalar_serialization.hpp>

// Make Mutation compatible with fwdpp's
// serialization API:
namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<fwdpy11::Mutation>
        {
            io::scalar_writer writer;
            serialize_mutation<fwdpy11::Mutation>() : writer{}
            {
            }
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const fwdpy11::Mutation &m) const
            {
                writer(buffer, &m.neutral);
                writer(buffer, &m.g);
                writer(buffer, &m.pos);
                writer(buffer, &m.s);
                writer(buffer, &m.h);
                writer(buffer, &m.xtra);
                std::size_t ns = m.esizes.size(), nh = m.heffects.size();
                writer(buffer, &ns);
                writer(buffer, &nh);
                if (ns)
                    {
                        writer(buffer, m.esizes.data(), ns);
                    }
                if (nh)
                    {
                        writer(buffer, m.heffects.data(), nh);
                    }
            }
        };

        template <> struct deserialize_mutation<fwdpy11::Mutation>
        {
            io::scalar_reader reader;
            deserialize_mutation<fwdpy11::Mutation>() : reader{}
            {
            }
            template <typename streamtype>
            inline fwdpy11::Mutation
            operator()(streamtype &buffer) const
            {
                std::int32_t g; // Changed from uint32_t in 0.13.0
                bool neutral;
                double pos, s, h;
                decltype(fwdpy11::Mutation::xtra) xtra;
                reader(buffer, &neutral);
                reader(buffer, &g);
                reader(buffer, &pos);
                reader(buffer, &s);
                reader(buffer, &h);
                reader(buffer, &xtra);
                std::size_t ns, nh;
                reader(buffer, &ns);
                reader(buffer, &nh);
                std::vector<double> ss, hs;
                if (ns)
                    {
                        ss.resize(ns);
                        reader(buffer, ss.data(), ns);
                    }
                if (nh)
                    {
                        hs.resize(ns);
                        reader(buffer, hs.data(), nh);
                    }
                return fwdpy11::Mutation(neutral, pos, s, h, g, std::move(ss),
                                         std::move(hs), xtra);
            }
        };
    }
}

#endif
