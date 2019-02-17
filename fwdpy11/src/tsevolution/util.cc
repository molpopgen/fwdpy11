#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/ts/definitions.hpp>

void
resize_genotype_matrix(std::vector<double> &new_diploid_gvalues,
                       std::size_t newsize, std::true_type)
{
    new_diploid_gvalues.resize(newsize);
}

void
resize_genotype_matrix(std::vector<double> & /*new_diploid_gvalues*/,
                       std::size_t /*newsize*/, std::false_type)
{
}

void
copy_genetic_values(double *beg, const std::vector<double> &gvalues,
                    std::true_type)
{
    std::copy(begin(gvalues), end(gvalues), beg);
}

void
copy_genetic_values(double * /*beg*/,
                    const std::vector<double> & /*genetic_value_fxn*/,
                    std::false_type)
{
}

void
remap_metadata(std::vector<fwdpy11::DiploidMetadata> &metadata,
               const std::vector<fwdpp::ts::TS_NODE_INT> &idmap)
{
    for (auto &m : metadata)
        {
            m.nodes[0] = idmap[m.nodes[0]];
            m.nodes[1] = idmap[m.nodes[1]];
            if (m.nodes[0] == fwdpp::ts::TS_NULL_NODE
                || m.nodes[1] == fwdpp::ts::TS_NULL_NODE)
                {
                    throw std::runtime_error(
                        "error remapping node field of individual metadata");
                }
        }
}

void
remap_ancient_samples(std::vector<fwdpy11::ancient_sample_record> &records,
                      const std::vector<fwdpp::ts::TS_NODE_INT> &idmap)
{
    for (auto &a : records)
        {
            a.n1 = idmap[a.n1];
            a.n2 = idmap[a.n2];
            if (a.n1 == fwdpp::ts::TS_NULL_NODE
                || a.n2 == fwdpp::ts::TS_NULL_NODE)
                {
                    throw std::runtime_error(
                        "error simplifying with respect to ancient samples");
                }
        }
}
