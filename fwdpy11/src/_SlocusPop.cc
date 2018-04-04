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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/SlocusPop.hpp>

namespace py = pybind11;

namespace
{
    static const auto DIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.VecDiploid`.
   )delim";
}

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Diploid>);

PYBIND11_MODULE(_SlocusPop, m)
{
    py::object base_class_module
        = (pybind11::object)pybind11::module::import("fwdpy11._Population");

    py::class_<fwdpy11::SlocusPop, fwdpy11::Population>(m, "_SlocusPop")
        .def(py::init<fwdpp::uint_t>())
        .def_readonly("diploids", &fwdpy11::SlocusPop::diploids,
                      DIPLOIDS_DOCSTRING);
}
