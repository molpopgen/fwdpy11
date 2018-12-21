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

#ifndef FWDPY11_GENETIC_VALUES_DETAILS_PICKLE_ADDITIVE_HPP
#define FWDPY11_GENETIC_VALUES_DETAILS_PICKLE_ADDITIVE_HPP

#include <pybind11/pybind11.h>
#include <fwdpp/fitness_models.hpp>

namespace fwdpy11
{
    struct pickle_additive
    {
        inline pybind11::object
        operator()(const fwdpp::additive_diploid& a) const
        {
            int p = 0;
            if (a.gvalue_is_trait == true)
                {
                    p = 1;
                }
            return pybind11::make_tuple(p, a.scaling);
        }
    };
} // namespace fwdpy11

#endif
