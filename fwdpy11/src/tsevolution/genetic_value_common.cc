#include <vector>
#include <cstdint>

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

