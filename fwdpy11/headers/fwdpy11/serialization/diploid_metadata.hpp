#ifndef FWDPY11_SERIALIZATION_DIPLOID_METADATA_HPP__
#define FWDPY11_SERIALIZATION_DIPLOID_METADATA_HPP__

#include <fwdpy11/types/Diploid.hpp>
#include <pybind11/pybind11.h>

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
            std::size_t s = vmd.size();
            w(buffer, &s);
            for (const auto& md : vmd)
                // For some reaason, serializing the
                // sub-arrays was causing segfaults.
                // TODO: revisit this issue later.
                {
                    w(buffer, &md.g);
                    w(buffer, &md.e);
                    w(buffer, &md.w);
                    w(buffer, &md.label);
                    w(buffer, &md.parents[0]);
                    w(buffer, &md.parents[1]);
                    w(buffer, &md.deme);
                    w(buffer, &md.sex);
                    w(buffer, &md.nodes[0]);
                    w(buffer, &md.nodes[1]);
                }
            pybind11::print("done");
        }
    };

    struct deserialize_diploid_metadata
    {
        template <typename streamtype>
        inline void
        operator()(streamtype& buffer,
                   std::vector<fwdpy11::dip_metadata>& vmd) const
        {
            fwdpp::io::scalar_reader r;
            std::size_t s;
            r(buffer, &s);
            pybind11::print(s);
            vmd.clear();

            dip_metadata md;
            for (std::size_t i = 0; i < s; ++i)
                {
                    r(buffer, &md.g);
                    r(buffer, &md.e);
                    r(buffer, &md.w);
                    r(buffer, &md.label);
                    r(buffer, &md.parents[0]);
                    r(buffer, &md.parents[1]);
                    r(buffer, &md.deme);
                    r(buffer, &md.sex);
                    r(buffer, &md.nodes[0]);
                    r(buffer, &md.nodes[1]);
                    vmd.emplace_back(md);
                }
        }
    };
} // namespace fwdpy11

#endif
