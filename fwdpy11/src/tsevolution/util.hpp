#ifndef FWDPY11_TSEVOLVE_UTIL_HPP
#define FWDPY11_TSEVOLVE_UTIL_HPP

#include <cstdint>
#include <vector>
#include <fwdpy11/types/Diploid.hpp>
#include <fwdpp/ts/definitions.hpp>

void resize_genotype_matrix(std::vector<double> &new_diploid_gvalues,
                            std::size_t newsize, std::true_type);

void resize_genotype_matrix(std::vector<double> & /*new_diploid_gvalues*/,
                            std::size_t /*newsize*/, std::false_type);

void copy_genetic_values(double *beg, const std::vector<double> &gvalues,
                         std::true_type);

void copy_genetic_values(double * /*beg*/,
                         const std::vector<double> & /*genetic_value_fxn*/,
                         std::false_type);

void
remap_metadata(std::vector<fwdpy11::DiploidMetadata> &metadata,
               const std::vector<fwdpp::ts::TS_NODE_INT> &idmap);

void
remap_ancient_samples(std::vector<fwdpy11::ancient_sample_record> &records,
                      const std::vector<fwdpp::ts::TS_NODE_INT> &idmap);
#endif
