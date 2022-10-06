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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/forward_types.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::haploid_genome>);

void init_HaploidGenomeVector(py::module & m)
{
    py::bind_vector<std::vector<fwdpp::haploid_genome>>(
        m, "HaploidGenomeVector", py::module_local(false),
        "C++ representations of a list of "
        ":class:`fwdpy11.HaploidGenome`.  "
        "Typically, access will be read-only.")
        .def(py::pickle(
            [](const std::vector<fwdpp::haploid_genome>& haploid_genomes) {
                py::list rv;
                for (auto&& g : haploid_genomes)
                    {
                        rv.append(g);
                    }
                return rv;
            },
            [](const py::list l) {
                std::vector<fwdpp::haploid_genome> rv;
                for (auto&& li : l)
                    {
                        rv.push_back(li.cast<fwdpp::haploid_genome>());
                    }
                return rv;
            }));
}
