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

#include <string>
#include <stdexcept>
#include <numeric>
#include <fwdpp/forward_types_serialization.hpp>
#include <fwdpp/io/serialize_population.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpy11/serialization/diploid_metadata.hpp>
#include "serialization/backwards_compat.hpp"

namespace fwdpy11
{
    namespace serialization
    {
        inline constexpr int
        magic()
        {
            // Changed to 3 im 0.3.0
            // to handle genetic value matrices
            // Changed to 4 in 0.5.2 to explicitly
            // handle the mutation counts, so that we can dodge tree
            // sequence traversal.
            // Changed to 5 in 0.5.3 to handle pop.mcounts explicitly.
            // The reason is that the fwdpp back-end regenerates
            // them by traversing genomes, but we may have neutral
            // mutations not in the genomes.
            // Changed to 6 in 0.6.3 because we removed "ancient sample
            // records" that weren't being used and we changed the C++
            // constructor for Mutation.
            return 6;
        }

        template <typename streamtype, typename poptype>
        streamtype &
        serialize_details(streamtype &buffer, const poptype *pop)
        {
            buffer << "fp11";
            auto m = magic();
            buffer.write(reinterpret_cast<char *>(&m), sizeof(decltype(m)));
            buffer.write(reinterpret_cast<const char *>((&pop->generation)),
                         sizeof(unsigned));
            fwdpy11::serialize_diploid_metadata()(buffer, pop->diploid_metadata);
            fwdpy11::serialize_diploid_metadata()(buffer, pop->ancient_sample_metadata);
            //fwdpy11::serialize_ancient_sample_records()(
            //   buffer, pop->ancient_sample_records);
            fwdpp::io::serialize_population(buffer, *pop);
            //preserved mutation counts added in 0.5.2, which is format version 4
            fwdpp::io::scalar_writer w;
            std::size_t msize = pop->mcounts_from_preserved_nodes.size();
            w(buffer, &msize);
            if (msize > 0)
                {
                    w(buffer, pop->mcounts_from_preserved_nodes.data(), msize);
                }
            // Added in 0.5.3/file version 5:
            msize = pop->mcounts.size();
            w(buffer, &msize);
            if (msize > 0)
                {
                    w(buffer, pop->mcounts.data(), msize);
                }
            fwdpp::ts::io::serialize_tables(buffer, *pop->tables);
            msize = pop->genetic_value_matrix.size();
            w(buffer, &msize);
            if (msize > 0)
                {
                    w(buffer, pop->genetic_value_matrix.data(), msize);
                }
            msize = pop->ancient_sample_genetic_value_matrix.size();
            w(buffer, &msize);
            if (msize > 0)
                {
                    w(buffer, pop->ancient_sample_genetic_value_matrix.data(), msize);
                }

            return buffer;
        }

        struct deserialize_details
        {
            template <typename streamtype, typename poptype>
            inline streamtype &
            operator()(streamtype &buffer, poptype &pop)
            {
                char c[4];
                buffer.read(&c[0], 4 * sizeof(char));
                std::string s(c, c + 4);
                // We need to test for existance of serialization
                // version numbers, introduced in 0.1.3.  Prior to that,
                // it was wild west :).
                bool have_magic = (s.substr(0, 4) == "fp11") ? true : false;
                // We default to version 1, which includes all
                // previous releases that had no version numbers
                int version = 1;
                if (have_magic)
                    {
                        buffer.read(reinterpret_cast<char *>(&version), sizeof(int));
                    }
                if (version == 1)
                    {
                        throw std::runtime_error("File format "
                                                 "incompatibility: this file "
                                                 "format version "
                                                 "was last supported in "
                                                 "fwdpy11 0.1.4");
                    }
                buffer.read(reinterpret_cast<char *>(&pop.generation), sizeof(unsigned));
                deserialize_diploid_metadata()(buffer, pop.diploid_metadata);
                deserialize_diploid_metadata()(buffer, pop.ancient_sample_metadata);
                if (version < 6)
                    {
                        std::vector<
                            deserialize_ancient_sample_records::ancient_sample_record>
                            ancient_sample_records;
                        fwdpy11::deserialize_ancient_sample_records()(
                            buffer, ancient_sample_records);
                        backwards_compat::deserialize_population_details(pop, buffer);
                    }
                else
                    {
                        fwdpp::io::deserialize_population(buffer, pop);
                    }
                std::size_t msize;
                fwdpp::io::scalar_reader r;
                if (version >= 4) // >= version 0.5.2
                    {
                        r(buffer, &msize);
                        pop.mcounts_from_preserved_nodes.resize(msize);
                        if (msize > 0)
                            {
                                r(buffer, pop.mcounts_from_preserved_nodes.data(),
                                  msize);
                            }
                    }
                if (version >= 5) // 0.5.3
                    {
                        r(buffer, &msize);
                        pop.mcounts.resize(msize);
                        if (msize > 0)
                            {
                                r(buffer, pop.mcounts.data(), msize);
                            }
                    }

                auto _tables = fwdpp::ts::io::deserialize_tables<fwdpp::ts::std_table_collection>()(buffer);
                pop.tables.reset(new fwdpp::ts::std_table_collection(std::move(_tables)));
                // NOTE: version 0.5.0 added in a site table that previous versions
                // did not have.  Further, the mutation table entries differed in
                // previous versions.
                if (!pop.tables->mutations.empty() && pop.tables->sites.empty())
                    {
                        fwdpp::ts::io::fix_mutation_table_repopulate_site_table(
                            *pop.tables, pop.mutations);
                    }
                if (version < 4 && !pop.tables->edges.empty())
                    {
                        std::vector<fwdpp::ts::table_index_t> samples(2 * pop.N);
                        std::iota(samples.begin(), samples.end(), 0);
                        fwdpp::ts::count_mutations(*pop.tables, pop.mutations, samples,
                                                   pop.mcounts);
                        pop.mcounts_from_preserved_nodes.resize(pop.mutations.size());
                        std::fill(begin(pop.mcounts_from_preserved_nodes),
                                  end(pop.mcounts_from_preserved_nodes), 0);
                    }
                if (version > 2)
                    {
                        r(buffer, &msize);
                        if (msize > 0)
                            {
                                pop.genetic_value_matrix.resize(msize);
                                r(buffer, pop.genetic_value_matrix.data(), msize);
                            }
                        r(buffer, &msize);
                        if (msize > 0)
                            {
                                pop.ancient_sample_genetic_value_matrix.resize(msize);
                                r(buffer, pop.ancient_sample_genetic_value_matrix.data(),
                                  msize);
                            }
                    }
                pop.rebuild_mutation_lookup(false);
                return buffer;
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
    } // namespace serialization
} // namespace fwdpy11
#endif
