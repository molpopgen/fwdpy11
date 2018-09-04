// clang-format off
<% 
setup_pybind11(cfg) 
cfg['libraries'].extend(['gsl','gslcblas'])
%>
// clang-format on

#include <pybind11/pybind11.h>
#include <gsl/gsl_matrix.h>

    namespace py = pybind11;

PYBIND11_MODULE(gsl_error, m)
{
    m.def("trigger_error", []() {
        gsl_matrix* m = gsl_matrix_alloc(2, 3); //non-square
        int rv = gsl_matrix_transpose(m);
        if (rv != GSL_SUCCESS)
            {
                throw std::runtime_error("non-square matrix");
                gsl_matrix_free(m);
            }
        gsl_matrix_free(m);
    });

    // Causes a memory leak!
    m.def("trigger_error_not_handled", []() {
        gsl_matrix* m = gsl_matrix_alloc(2, 3); //non-square
        gsl_matrix_transpose(m);
        gsl_matrix_free(m);
    });
}
