#ifndef FWDPY11_GSL_GSL_ERROR_HANDLER_WRAPPER_HPP
#define FWDPY11_GSL_GSL_ERROR_HANDLER_WRAPPER_HPP

#include <gsl/gsl_errno.h>

namespace fwdpy11
{
    struct gsl_error_handler_wrapper
    {
        gsl_error_handler_t* e;
        gsl_error_handler_wrapper() : e{ gsl_set_error_handler_off() } {}
        ~gsl_error_handler_wrapper() { gsl_set_error_handler(e); }
    };
} // namespace fwdpy11

#endif
