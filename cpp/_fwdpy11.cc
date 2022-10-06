#include <stdexcept>
#include <type_traits>
#include <gsl/gsl_version.h>
#include <gsl/gsl_errno.h>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
#include <pybind11/pybind11.h>

static_assert(GSL_MAJOR_VERSION >= 2, "GSL major version >= 2 required");
static_assert(GSL_MINOR_VERSION >= 3, "GSL minor version >= 3 required");

namespace py = pybind11;

void initialize_fwdpp_types(py::module &);
void initialize_fwdpy11_types(py::module &m);
void initialize_regions(py::module &);
void initialize_mutation_dominance(py::module &m);
void initialize_genetic_value_noise(py::module &);
void initialize_genetic_value_to_fitness(py::module &);
void init_genetic_values(py::module &);
void init_GSL(py::module &);
void init_ts(py::module &);
void init_evolution_functions(py::module &);
void init_discrete_demography(py::module &m);
void init_array_proxies(py::module &m);
void initialize_functions(py::module &m);
void init_demes(py::module &m);

PYBIND11_MODULE(_fwdpy11, m)
{
    initialize_fwdpp_types(m);
    initialize_fwdpy11_types(m);
    initialize_mutation_dominance(m);
    initialize_regions(m);
    initialize_genetic_value_noise(m);
    initialize_genetic_value_to_fitness(m);
    init_genetic_values(m);
    init_GSL(m);
    init_ts(m);
    init_evolution_functions(m);
    init_discrete_demography(m);
    init_array_proxies(m);
    initialize_functions(m);
    init_demes(m);

    py::register_exception<fwdpy11::GSLError>(m, "GSLError");

    m.def(
        "pybind11_version",
        []() {
            py::dict rv;
            std::ostringstream o;
            o << PYBIND11_VERSION;
            rv["pybind11_version"] = o.str();
            return rv;
        },
        R"delim(
    Returns the version of pybind11 used to
    compile fwdpy11.
    )delim");
}
