#ifndef FWDPY11_SAMPLING_DATA_MATRIX_FUNCTIONS_HPP__
#define FWDPY11_SAMPLING_DATA_MATRIX_FUNCTIONS_HPP__

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <array>
#include <utility>
#include <string>
#include <fwdpp/data_matrix.hpp>
#include <gsl/gsl_matrix_char.h>

namespace fwdpy11
{
    std::vector<std::pair<double, std::string>>
    matrix_to_sample(const fwdpp::state_matrix &sm)
    // returns a data structure compatible with libsequence/pylibseq iff
    // the data correspond to a haplotype matrix
    {
        std::vector<std::pair<double, std::string>> rv;
        if (sm.positions.empty())
            {
                return rv;
            }
        if (sm.data.empty())
            {
                throw std::runtime_error("StateMatrix data are empty");
            }
        std::size_t nrow = sm.positions.size();
        std::size_t ncol = sm.data.size()/nrow;
        const std::array<std::int8_t, 3> states{ '0', '1', '2' };
        auto v = gsl_matrix_char_const_view_array(
            reinterpret_cast<const char *>(sm.data.data()), nrow, ncol);
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
                rv.push_back(std::make_pair(sm.positions[i], std::move(column_data)));
            }
        return rv;
    }
} // namespace fwdpy11
#endif
