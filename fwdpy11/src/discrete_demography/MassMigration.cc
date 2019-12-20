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
#include <sstream>
#include <pybind11/pybind11.h>
#include <fwdpy11/discrete_demography/MassMigration.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

ddemog::MassMigration
move_individuals(std::uint32_t when, std::int32_t source,
                 std::int32_t destination, double fraction,
                 bool resets_growth_rate)
{
    return ddemog::MassMigration(source, destination, when, 0, -1, fraction,
                                  true, false, resets_growth_rate);
}

ddemog::MassMigration
copy_individuals(std::uint32_t when, std::int32_t source,
                 std::int32_t destination, double fraction,
                 bool resets_growth_rate)
{
    return ddemog::MassMigration(source, destination, when, 0, -1, fraction,
                                  false, false, resets_growth_rate);
}

static const auto MOVE_INDIVIDUALS_DOCSTRING = R"delim(
:param when: The generation when the event occurs
:type when: int
:param source: The source deme for individuals to move
:type source: int
:param destination: The deme to where individuals we be moved
:type destination: int
:param fraction: The fraction of `source` to move to `destination`.
:type fraction: float
:param resets_growth_rate: (True) Whether or not to reset the growth rates of `source` and `destination` to :data:`fwdpy11.NOGROWTH`
:type resets_growth_rate: bool

:rtype: :class:`fwdpy11.MassMigration`
)delim";

static const auto COPY_INDIVIDUALS_DOCSTRING = R"delim(
:param when: The generation when the event occurs
:type when: int
:param source: The source deme for individuals to copy
:type source: int
:param destination: The deme to where individuals we be copied
:type destination: int
:param fraction: The fraction of `source` to copy to `destination`.
:type fraction: float
:param resets_growth_rate: (True) Whether or not to reset the growth rates of `source` and `destination` to :data:`fwdpy11.NOGROWTH`
:type resets_growth_rate: bool

:rtype: :class:`fwdpy11.MassMigration`
)delim";

namespace
{
    std::string
    repr(const ddemog::MassMigration& self)
    {
        std::ostringstream o;
        o << "MassMigration(when=" << self.when << ", source=" << self.source
          << ", destination=" << self.destination
          << ", fraction=" << self.fraction
          << ", move_individuals=" << self.move_individuals << ')';
        return o.str();
    }
} // namespace

void
init_MassMigration(py::module& m)
{
    py::class_<ddemog::MassMigration>(m, "MassMigration",
                                       R"delim(
        Mass migration events.

        .. versionadded:: 0.5.3
        )delim")
        .def_readonly("when", &ddemog::MassMigration::when)
        .def_readonly("source", &ddemog::MassMigration::source)
        .def_readonly("destination", &ddemog::MassMigration::destination)
        .def_readonly("number", &ddemog::MassMigration::number)
        .def_readonly("fraction", &ddemog::MassMigration::fraction)
        .def_readonly("move_individuals",
                      &ddemog::MassMigration::move_individuals)
        .def_readonly("resets_growth_rate",
                      &ddemog::MassMigration::resets_growth_rate)
        .def("__repr__",
             [](const ddemog::MassMigration& self) { return repr(self); })
        .def(py::pickle(
            [](const ddemog::MassMigration& self) {
                return py::make_tuple(self.source, self.destination, self.when,
                                      self.number, self.sex, self.fraction,
                                      self.move_individuals, self.sex_specific,
                                      self.resets_growth_rate);
            },
            [](py::tuple t) {
                using M = ddemog::MassMigration;
                return ddemog::MassMigration(
                    t[0].cast<decltype(M::source)>(),
                    t[1].cast<decltype(M::destination)>(),
                    t[2].cast<decltype(M::when)>(),
                    t[3].cast<decltype(M::number)>(),
                    t[4].cast<decltype(M::sex)>(),
                    t[5].cast<decltype(M::fraction)>(),
                    t[6].cast<decltype(M::move_individuals)>(),
                    t[7].cast<decltype(M::sex_specific)>(),
                    t[8].cast<decltype(M::resets_growth_rate)>());
            }))
        .def("__eq__",
             [](const ddemog::MassMigration& lhs,
                const ddemog::MassMigration& rhs) { return lhs == rhs; });

    m.def("move_individuals", &move_individuals, py::arg("when"),
          py::arg("source"), py::arg("destination"), py::arg("fraction"),
          py::arg("resets_growth_rate") = true, MOVE_INDIVIDUALS_DOCSTRING);
    m.def("copy_individuals", &copy_individuals, py::arg("when"),
          py::arg("source"), py::arg("destination"), py::arg("fraction"),
          py::arg("resets_growth_rate") = true, COPY_INDIVIDUALS_DOCSTRING);
}
