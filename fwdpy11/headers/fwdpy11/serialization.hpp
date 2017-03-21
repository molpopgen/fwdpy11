/*! \file fwdpy_serialization.hpp
 * \brief Helper functions for object-level serialization
 */
#ifndef FWDPY_SERIALIZATION_HPP
#define FWDPY_SERIALIZATION_HPP
#include "serialization_common.hpp"
#include <fwdpp/sugar/serialization.hpp>
namespace fwdpy
{
    namespace serialize_objects
    {

        template <typename poptype> struct deserialize_details
        {
            template <typename mreader_t, typename dipreader_t,
                      typename... constructor_data>
            inline poptype
            operator()(const std::string &s, const mreader_t &mreader,
                       const dipreader_t &dipreader, constructor_data... cdata)
            {
                std::istringstream buffer;
                buffer.str(s);
                buffer.seekg(0);
                poptype pop(cdata...);
                buffer.read(reinterpret_cast<char *>(&pop.generation),
                            sizeof(unsigned));
                KTfwd::deserialize d;
                d(pop, buffer, mreader, dipreader);
                return pop;
            }
        };

        template <typename poptype, typename mwriter_t, typename dipwriter_t>
        inline int
        gzserialize_details(const poptype &pop, const mwriter_t &mwriter,
                            const dipwriter_t &dipwriter, const char *filename,
                            bool append)
        {
            gzFile f;
            if (append)
                {
                    f = gzopen(filename, "ab");
                }
            else
                {
                    f = gzopen(filename, "wb");
                }
            auto rv
                = gzwrite(f, reinterpret_cast<const char *>(&pop.generation),
                          sizeof(decltype(pop.generation)));
            KTfwd::gzserialize s;
            rv += s(f, pop, mwriter, dipwriter);
            gzclose(f);
            return rv;
        }

        template <typename poptype> struct gzdeserialize_details
        {
            template <typename mreader_t, typename dipreader_t,
                      typename... constructor_data>
            inline poptype
            operator()(const mreader_t &mreader, const dipreader_t &dipreader,
                       const char *filename, std::size_t offset,
                       constructor_data... cdata) const
            {
                gzFile f = gzopen(filename, "rb");
                if (offset)
                    {
                        gzseek(f, offset, SEEK_SET);
                    }
                poptype temp(cdata...);
                gzread(f, reinterpret_cast<char *>(&temp.generation),
                       sizeof(decltype(temp.generation)));
                KTfwd::gzdeserialize s;
                s(temp, f, mreader, dipreader);
                gzclose(f);
                return temp;
            };
        };
    }
}
#endif
