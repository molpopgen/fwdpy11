#ifndef FWDPY_SERIALIATION_COMMON_HPP
#define FWDPY_SERIALIATION_COMMON_HPP
#include <fwdpp/sugar/serialization.hpp>
#include <sstream>
#include <string>
namespace fwdpy
{
    namespace serialization
    {

        template <typename poptype, typename mwriter_t, typename dipwriter_t>
        std::string
        serialize_details(const poptype *pop, const mwriter_t &mwriter,
                          const dipwriter_t &dipwriter)
        {
            KTfwd::serialize rv;
            std::ostringstream buffer;
            buffer.write(reinterpret_cast<const char *>((&pop->generation)),
                         sizeof(unsigned));
            rv(buffer, *pop, mwriter, dipwriter);
            return buffer.str();
        }
    }
}

#endif
