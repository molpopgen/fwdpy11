#ifndef FWDPY11_SERIALIZATION_DIPLOID_METADATA_HPP__
#define FWDPY11_SERIALIZATION_DIPLOID_METADATA_HPP__

#include <fwdpy11/types/Diploid.hpp>

namespace fwdpy11
{
    struct serialize_diploid_metadata
    {
        template <typename streamtype>
        inline void
        operator()(streamtype& buffer,
                   const std::vector<fwdpy11::dip_metadata>& vmd) const
        {
            fwdpp::io::scalar_writer w;
            auto s = vmd.size();
            w(buffer, &s);
            for (const auto& md : vmd)
                {
                    w(buffer, &md.g);
                    w(buffer, &md.e);
                    w(buffer, &md.w);
                    w(buffer, &md.parents, 2);
                    w(buffer, &md.deme);
                    w(buffer, &md.sex);
                    w(buffer, &md.nodes, 2);
                }
        }
    };

    struct deserialize_diploid_metadata
    {
        template <typename streamtype>
        inline void
        operator()(streamtype& buffer,
                   std::vector<fwdpy11::dip_metadata>& vmd) const
        {
            std::size_t s;
            fwdpp::io::scalar_reader r;
            r(buffer, &s);
            vmd.clear();
            fwdpy11::dip_metadata md;
            for (std::size_t i = 0; i < s; ++i)
                {
                    r(buffer, &md.g);
                    r(buffer, &md.e);
                    r(buffer, &md.w);
                    r(buffer, &md.parents, 2);
                    r(buffer, &md.deme);
                    r(buffer, &md.sex);
                    r(buffer, &md.nodes, 2);
                    vmd.push_back(md);
                }
        }
    };
} // namespace fwdpy11


#endif
