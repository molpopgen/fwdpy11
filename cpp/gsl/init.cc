#include <sstream>
#include <gsl/gsl_version.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_gsl_random(py::module&);

void
init_GSL(py::module& m)
{
    init_gsl_random(m);

    m.def(
        "gsl_version",
        []() {
            py::dict rv;
            std::ostringstream o;
            o << GSL_MAJOR_VERSION << '.' << GSL_MINOR_VERSION;
            rv["gsl_version"] = o.str();
            return rv;
        },
        R"delim(
    Returns the version of the GSL used to compile fwdpy11.

    :rtype: dict
    
    .. versionadded:: 0.5.0
    )delim");
}
