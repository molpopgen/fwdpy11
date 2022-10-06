//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

#include <fwdpy11/util/array_proxy.hpp>

namespace py = pybind11;

void
init_array_proxies(py::module& m)
{
    py::class_<fwdpy11::double_array_proxy>(m, "_FloatProxy", py::buffer_protocol())
        .def_buffer(
            [](const fwdpy11::double_array_proxy& self) { return as_buffer(self); });

    py::class_<fwdpy11::uint32_array_proxy>(m, "_Uint32ArrayProxy",
                                            py::buffer_protocol())
        .def_buffer(
            [](const fwdpy11::uint32_array_proxy& self) { return as_buffer(self); });
}
