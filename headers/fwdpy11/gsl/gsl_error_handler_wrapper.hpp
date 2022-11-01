#ifndef FWDPY11_GSL_GSL_ERROR_HANDLER_WRAPPER_HPP
#define FWDPY11_GSL_GSL_ERROR_HANDLER_WRAPPER_HPP

#include <exception>
#include <sstream>
#include <string>
#include <gsl/gsl_errno.h>

namespace fwdpy11
{

    class __attribute__((visibility("default"))) GSLError : public std::exception
    {
      private:
        std::string message_;

      public:
        explicit GSLError(std::string message) : message_(std::move(message))
        {
        }
        virtual const char*
        what() const noexcept
        {
            return message_.c_str();
        }
    };

    class gsl_scoped_convert_error_to_exception
    /*!
     * Manages turning off and resetting the GSL error
     * handler.
     *
     * The handler is restored via the destructor,
     * so this class is like a "smart pointer" for
     * turning off the handler.
     */
    {
      private:
        gsl_error_handler_t *e, *custom_handler;

        static void
        gsl_error_to_exception(const char* reason, const char* file, int line,
                               int gsl_errno)
        {
            std::ostringstream o;
            o << "GSL error raised: " << reason << ", " << file << ", " << line << ", "
              << gsl_errno;
            throw fwdpy11::GSLError(o.str());
        }

      public:
        gsl_scoped_convert_error_to_exception()
            : e{gsl_set_error_handler_off()},
              custom_handler{gsl_set_error_handler(
                  &gsl_scoped_convert_error_to_exception::gsl_error_to_exception)}
        {
        }

        ~gsl_scoped_convert_error_to_exception()
        {
            gsl_set_error_handler(e);
        }

        gsl_scoped_convert_error_to_exception(
            const gsl_scoped_convert_error_to_exception&)
            = delete;
        gsl_scoped_convert_error_to_exception(gsl_scoped_convert_error_to_exception&&)
            = default;
        gsl_scoped_convert_error_to_exception&
        operator=(const gsl_scoped_convert_error_to_exception&)
            = delete;
        gsl_scoped_convert_error_to_exception&
        operator=(gsl_scoped_convert_error_to_exception&&)
            = default;
    };

    struct gsl_scoped_disable_error_handler_wrapper
    {
        gsl_error_handler_t* e;

        gsl_scoped_disable_error_handler_wrapper() : e{gsl_set_error_handler_off()}
        {
        }

        ~gsl_scoped_disable_error_handler_wrapper()
        {
            gsl_set_error_handler(e);
        }

        gsl_scoped_disable_error_handler_wrapper(
            const gsl_scoped_disable_error_handler_wrapper&)
            = delete;
        gsl_scoped_disable_error_handler_wrapper(
            gsl_scoped_disable_error_handler_wrapper&&)
            = default;

        gsl_scoped_disable_error_handler_wrapper&
        operator=(const gsl_scoped_disable_error_handler_wrapper&)
            = delete;
        gsl_scoped_disable_error_handler_wrapper&
        operator=(gsl_scoped_disable_error_handler_wrapper&&)
            = default;
    };
} // namespace fwdpy11

#endif
