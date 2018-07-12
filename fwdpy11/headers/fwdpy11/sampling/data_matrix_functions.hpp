#ifndef FWDPY11_SAMPLING_DATA_MATRIX_FUNCTIONS_HPP__
#define FWDPY11_SAMPLING_DATA_MATRIX_FUNCTIONS_HPP__

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <array>
#include <utility>
#include <string>
#include <gsl/gsl_matrix_char.h>

namespace fwdpy11
{
    std::vector<std::pair<double, std::string>>
    matrix_to_sample(const std::vector<std::int8_t> &data,
                     const std::vector<double> &pos, const std::size_t ncol)
    // returns a data structure compatible with libsequence/pylibseq iff
    // the data correspond to a haplotype matrix
    {
        std::size_t nrow = data.size() / ncol;
        const std::array<std::int8_t, 3> states{ '0', '1', '2' };
        auto v = gsl_matrix_char_const_view_array(
            reinterpret_cast<const char *>(data.data()), nrow, ncol);
        std::vector<std::pair<double, std::string>> rv;
        for (std::size_t i = 0; i < nrow; ++i)
            {
                auto c = gsl_matrix_char_const_row(&v.matrix, i);
                std::string column_data;
                for (std::size_t j = 0; j < c.vector.size; ++j)
                    {
                        column_data.push_back(states[static_cast<std::int8_t>(
                            gsl_vector_char_get(&c.vector, j))]);
                    }
                if (column_data.size() != ncol)
                    {
                        throw std::runtime_error("column_data.size() != ncol");
                    }
                rv.push_back(std::make_pair(pos[i], std::move(column_data)));
            }
        return rv;
    }

    std::vector<std::vector<std::pair<double, std::string>>>
    separate_samples_by_loci(
        const std::vector<std::pair<double, double>> &boundaries,
        const std::vector<std::pair<double, std::string>> &sample)
    // For a multi-locus pop, it is convenient to split samples by
    // loci.  This function does that using pop.locus_boundaries.
    // If pop.locus_boundaries is not properly set, an exception
    // is likely going to be triggered
    // The parameter "sample" is the return value of matrix_to_sample.
    {
        std::vector<std::vector<std::pair<double, std::string>>> rv(boundaries.size());
        if (sample.size() == 0)
            {
                return rv;
            }
        for (auto &&site : sample)
            {
                auto itr = std::find_if(
                    boundaries.begin(), boundaries.end(),
                    [&site](const std::pair<double, double> &b) {
                        return site.first >= b.first
                               && site.first < b.second;
                    });
                if (itr == boundaries.end())
                    {
                        throw std::runtime_error(
                            "could not find locus for mutation at position"
                            + std::to_string(site.first));
                    }
                auto d = std::distance(boundaries.begin(), itr);
                rv[d].push_back(site);
            }
        return rv;
    }
} // namespace fwdpy11
#endif
