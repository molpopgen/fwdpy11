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
#include "DiploidGeneticValue.hpp"

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
        std::vector<fwdpy11::DiploidGeneticValue *> genetic_values;

        dgvalue_pointer_vector_(fwdpy11::DiploidGeneticValue &gv) : genetic_values{&gv}
        {
        }

        dgvalue_pointer_vector_(std::vector<fwdpy11::DiploidGeneticValue *> gvpointers)
            : genetic_values{std::move(gvpointers)}
        {
        }
    };
} // namespace fwdpy11

#endif
