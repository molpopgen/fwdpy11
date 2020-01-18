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

#ifndef FWDPY11_GENETIC_VALUES_DGVALUE_POINTER_VECTOR
#define FWDPY11_GENETIC_VALUES_DGVALUE_POINTER_VECTOR

#include <vector>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include "DiploidPopulationGeneticValue.hpp"

namespace fwdpy11
{
    struct dgvalue_pointer_vector_
    // For internal use only.
    // Introduced in 0.6.0 to manage
    // vectors of genetic value objecst in multi-deme
    // simulations.
    // This class exists because our genetic value types
    // are held by Python in std::unique_ptr, and
    // pybind11 (correctly) doesn't allow unique_ptr<T>
    // to be passed as an argument to a C++ function.
    // Thus, we compromise with raw pointers stored in a vector.
    {
        const std::vector<fwdpy11::DiploidPopulationGeneticValue *>
            genetic_values;

        std::vector<fwdpy11::DiploidPopulationGeneticValue *>
        init_from_list(pybind11::list l)
        {
            if (l.empty())
                {
                    throw std::invalid_argument(
                        "list of genetic values cannot be empty");
                }
            std::vector<fwdpy11::DiploidPopulationGeneticValue *> rv;
            for (auto i : l)
                {
                    auto *ref
                        = i.cast<fwdpy11::DiploidPopulationGeneticValue *>();
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

        dgvalue_pointer_vector_(fwdpy11::DiploidPopulationGeneticValue &gv)
            : genetic_values{ &gv }
        {
        }

        dgvalue_pointer_vector_(pybind11::list gv)
            : genetic_values{ init_from_list(gv) }
        {
        }
    };
} // namespace fwdpy11

#endif
