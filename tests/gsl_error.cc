#include <pybind11/pybind11.h>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
#include <gsl/gsl_matrix.h>

PYBIND11_MODULE(gsl_error, m)
{
    m.def("trigger_error_handle_locally", []() {
        fwdpy11::gsl_scoped_disable_error_handler_wrapper h;
        gsl_matrix* m = gsl_matrix_alloc(2, 3); //non-square
        int rv = gsl_matrix_transpose(m);
        if (rv != GSL_SUCCESS)
            {
                gsl_matrix_free(m);
                throw std::runtime_error("matrix not square");
            }
        gsl_matrix_free(m);
    });

    m.def("trigger_error", []() {
        fwdpy11::gsl_scoped_convert_error_to_exception gsl_error_scope;
        gsl_matrix* m = gsl_matrix_alloc(2, 3); //non-square
        try
            {
                gsl_matrix_transpose(m);
                gsl_matrix_free(m);
            }
        catch (fwdpy11::GSLError& e)
            {
                throw e;
            }
        catch (...)
            {
                throw std::runtime_error("unexcpected exception");
            }
    });
}
