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
#include <type_traits>

static_assert(GSL_MAJOR_VERSION >= 2, "GSL major version >= 2 required");
static_assert(GSL_MINOR_VERSION >= 2, "GSL minor version >= 2 required");

PYBIND11_MODULE(_init, m)
{
    m.doc() = "Module executing some setup code when fwdpy11 is imported.";

    auto handler = gsl_set_error_handler_off();
}
