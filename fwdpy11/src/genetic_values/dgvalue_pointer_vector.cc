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
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>

namespace py = pybind11;

namespace
{
    std::vector<fwdpy11::DiploidGeneticValue *>
    init_from_list(pybind11::list l)
    {
        if (l.empty())
            {
                throw std::invalid_argument("list of genetic values cannot be empty");
            }
        std::vector<fwdpy11::DiploidGeneticValue *> rv;
        for (auto i : l)
            {
                auto *ref = i.cast<fwdpy11::DiploidGeneticValue *>();
                rv.push_back(ref);
            }
        for (std::size_t i = 1; i < rv.size(); ++i)
            {
                if (rv[i - 1]->total_dim != rv[i]->total_dim)
                    {
                        rv.clear();
                        throw std::invalid_argument(
                            "genetic value objects must all have same "
                            "value for total_dim");
                    }
            }
        return rv;
    }
}

void
init_dgvalue_pointer_vector(py::module &m)
{
    py::class_<fwdpy11::dgvalue_pointer_vector_>(m, "_dgvalue_pointer_vector")
        .def(py::init<fwdpy11::DiploidGeneticValue &>())
        .def(py::init([](py::list l) {
            auto pointers = init_from_list(l);
            return fwdpy11::dgvalue_pointer_vector_(std::move(pointers));
        }));
}
