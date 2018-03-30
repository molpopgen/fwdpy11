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
/*! \file fwdpy_serialization.hpp
 * \brief Helper functions for object-level serialization
 */
#ifndef FWDPY11_SERIALIZATION_HPP
#define FWDPY11_SERIALIZATION_HPP
#include <sstream>
#include <fwdpp/forward_types_serialization.hpp>
#include <fwdpp/io/serialize_population.hpp>

namespace fwdpy11
{
    namespace serialization
    {
        inline constexpr int
        magic()
        {
            return 2;
        }

        template <typename poptype>
        std::string
        serialize_details(const poptype *pop)
        {
            std::ostringstream buffer;
            buffer << "fp11";
            auto m = magic();
            buffer.write(reinterpret_cast<char *>(&m), sizeof(decltype(m)));
            buffer.write(reinterpret_cast<const char *>((&pop->generation)),
                         sizeof(unsigned));
            fwdpp::io::serialize_population(buffer, *pop);
            return buffer.str();
        }

        template <typename poptype> struct deserialize_details
        {
            template <typename... constructor_data>
            inline poptype
            operator()(const std::string &s, constructor_data... cdata)
            {
                std::istringstream buffer;
                buffer.str(s);
                buffer.seekg(0);
                poptype pop(cdata...);
                // We need to test for existance of serialization
                // version numbers, introduced in 0.1.3.  Prior to that,
                // it was wild west :).
                bool have_magic = (s.substr(0, 4) == "fp11") ? true : false;
                // We default to version 1, which includes all
                // previous releases that had no version numbers
                int version = 1;
                if (have_magic)
                    {
                        char c[4];
                        buffer.read(&c[0], 4 * sizeof(char));
                        buffer.read(reinterpret_cast<char *>(&version),
                                    sizeof(int));
                    }
                if (version == 1)
                    {
                        throw std::runtime_error("File format "
                                                 "incompatibility: this file "
                                                 "format version "
                                                 "was last supported in "
                                                 "fwdpy11 0.1.4");
                    }
                buffer.read(reinterpret_cast<char *>(&pop.generation),
                            sizeof(unsigned));
                fwdpp::io::deserialize_population(buffer, pop);
                return pop;
            }
        };

        // template <typename poptype, typename mwriter_t, typename
        // dipwriter_t> inline int gzserialize_details(const poptype &pop,
        // const mwriter_t &mwriter,
        //                     const dipwriter_t &dipwriter, const char
        //                     *filename, bool append)
        // {
        //     gzFile f;
        //     if (append)
        //         {
        //             f = gzopen(filename, "ab");
        //         }
        //     else
        //         {
        //             f = gzopen(filename, "wb");
        //         }
        //     auto rv
        //         = gzwrite(f, reinterpret_cast<const char
        //         *>(&pop.generation),
        //                   sizeof(decltype(pop.generation)));
        //     fwdpp::gzserialize s;
        //     rv += s(f, pop, mwriter, dipwriter);
        //     gzclose(f);
        //     return rv;
        // }

        // template <typename poptype> struct gzdeserialize_details
        // {
        //     template <typename mreader_t, typename dipreader_t,
        //               typename... constructor_data>
        //     inline poptype
        //     operator()(const mreader_t &mreader, const dipreader_t
        //     &dipreader,
        //                const char *filename, std::size_t offset,
        //                constructor_data... cdata) const
        //     {
        //         gzFile f = gzopen(filename, "rb");
        //         if (offset)
        //             {
        //                 gzseek(f, offset, SEEK_SET);
        //             }
        //         poptype temp(cdata...);
        //         gzread(f, reinterpret_cast<char *>(&temp.generation),
        //                sizeof(decltype(temp.generation)));
        //         fwdpp::gzdeserialize s;
        //         s(temp, f, mreader, dipreader);
        //         gzclose(f);
        //         return temp;
        //     };
        // };
    }
}
#endif
