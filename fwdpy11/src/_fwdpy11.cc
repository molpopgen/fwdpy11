#include <stdexcept>
#include <sstream>
#include <type_traits>
#include <gsl/gsl_version.h>
#include <gsl/gsl_errno.h>
#include <pybind11/pybind11.h>

static_assert(GSL_MAJOR_VERSION >= 2, "GSL major version >= 2 required");
static_assert(GSL_MINOR_VERSION >= 3, "GSL minor version >= 3 required");

void
gsl_error_to_exception(const char *reason, const char *file, int line,
                       int gsl_errno)
{
    std::ostringstream o;
    o << "GSL error raised: " << reason << ", " << file << ", " << line << ", "
      << gsl_errno;
    throw std::runtime_error(o.str());
}

namespace py = pybind11;

void initialize_fwdpp_types(py::module &);
void initialize_fwdpp_functions(py::module &);
void initialize_fwdpy11_types(py::module &m);
void initialize_fwdpy11_functions(py::module &);
void initialize_regions(py::module &);
void initialize_genetic_value_noise(py::module &);
void initialize_genetic_value_to_fitness(py::module &);
void init_genetic_values(py::module &);
void init_GSL(py::module &);
void init_ts(py::module&);
void init_evolution_functions(py::module&);
void init_discrete_demography(py::module & m);

PYBIND11_MODULE(_fwdpy11, m)
{
    auto handler = gsl_set_error_handler_off();

    auto custom_handler = gsl_set_error_handler(&gsl_error_to_exception);

    // When the module exits, restore the old handler.
    // We re-use variables to suppress unused variable warnings.
    py::module::import("atexit").attr("register")(
        py::cpp_function{ [&handler, &custom_handler]() {
            custom_handler = gsl_set_error_handler(handler);
        } });

    initialize_fwdpp_types(m);
    initialize_fwdpp_functions(m);
    initialize_fwdpy11_types(m);
    initialize_fwdpy11_functions(m);
    initialize_regions(m);
    initialize_genetic_value_noise(m);
    initialize_genetic_value_to_fitness(m);
    init_genetic_values(m);
    init_GSL(m);
    init_ts(m);
    init_evolution_functions(m);
    init_discrete_demography(m);
}
