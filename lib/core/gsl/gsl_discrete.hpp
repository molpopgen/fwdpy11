#pragma once

#include <fwdpp/gsl_discrete.hpp>

namespace fwdpy11_core
{
    void update_lookup_table(const double *data, std::size_t size,
                             fwdpp::gsl_ran_discrete_t_ptr &lookup);
    void update_lookup_table_skip_zero_check(const double *data, std::size_t size,
                                             fwdpp::gsl_ran_discrete_t_ptr &lookup);
}
