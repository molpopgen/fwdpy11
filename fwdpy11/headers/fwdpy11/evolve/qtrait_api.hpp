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
#ifndef FWDPY11_EVOLVE_QTRAIT_API__
#define FWDPY11_EVOLVE_QTRAIT_API__

#include <functional>
#include <fwdpy11/types.hpp>

namespace fwdpy11
{
    using trait_to_fitness_function
        = std::function<double(const double, const double)>;
    using single_locus_noise_function = std::function<double(
        const fwdpy11::diploid_t &, const fwdpy11::diploid_t &)>;
    using multilocus_noise_function = std::function<double(
        const fwdpy11::multilocus_diploid_t &, const fwdpy11::multilocus_diploid_t &)>;
    using multilocus_aggregator_function
        = std::function<double(const pybind11::array_t<double>)>;
}

#endif
