#include <algorithm>

#include <core/gsl/gsl_discrete.hpp>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>

namespace fwdpy11_core
{
    void
    update_lookup_table(const double *data, std::size_t size,
                        fwdpp::gsl_ran_discrete_t_ptr &lookup)
    {
        if (std::all_of(data, data + size, [](const double d) { return d == 0.0; }))
            {
                throw fwdpy11::GSLError("all weights are 0.0");
            }
        update_lookup_table_skip_zero_check(data, size, lookup);
    }

    void
    update_lookup_table_skip_zero_check(const double *data, std::size_t size,
                                        fwdpp::gsl_ran_discrete_t_ptr &lookup)
    {
        // Note, without a handler from
        // #include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
        // the following call with abort when presented invalid data
        auto table = gsl_ran_discrete_preproc(size, data);
        if (table == nullptr)
            {
                throw fwdpy11::GSLError("gsl_ran_discrete_preproc returned nullptr");
            }
        lookup.reset(table);
    }
}
