//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <pybind11/pybind11.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <sstream>
#include <type_traits>

static_assert(GSL_MAJOR_VERSION >= 2, "GSL major version >= 2 required");
static_assert(GSL_MINOR_VERSION >= 3, "GSL minor version >= 3 required");

namespace py = pybind11;

void
gsl_error_to_exception(const char* reason, const char* file, int line,
                       int gsl_errno)
{
    std::ostringstream o;
    o << "GSL error raised: " << reason << ", " << file << ", " << line << ", "
      << gsl_errno;
    throw std::runtime_error(o.str());
}

PYBIND11_MODULE(_init, m)
{
    m.doc() = "Module executing some setup code when fwdpy11 is imported.";

    auto handler = gsl_set_error_handler_off();

    auto custom_handler = gsl_set_error_handler(&gsl_error_to_exception);

    // When the module exits, restore the old handler.
    // We re-use variables to suppress unused variable warnings.
    py::module::import("atexit").attr("register")(
        py::cpp_function{ [&handler, &custom_handler]() {
            custom_handler = gsl_set_error_handler(handler);
        } });
}
