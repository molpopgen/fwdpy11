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
#ifndef FWDPY11_GENETIC_VALUES_MULT_HPP__
#define FWDPY11_GENETIC_VALUES_MULT_HPP__

#include <fwdpp/fitness_models.hpp>
#include "details/pickle_multiplicative.hpp"
#include "fwdpp_wrappers/fwdpp_genetic_value.hpp"

namespace fwdpy11
{
    using DiploidMult = fwdpp_genetic_value<fwdpp::multiplicative_diploid,
                                           pickle_multiplicative>;
} // namespace fwdpy11
#endif
