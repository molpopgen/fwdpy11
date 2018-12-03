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
#ifndef FWDPY11_GENETIC_VALUES_MLOCUSADDITIVE_HPP__
#define FWDPY11_GENETIC_VALUES_MLOCUSADDITIVE_HPP__

#include "fwdpp_wrappers/fwdpp_mlocus_gvalue.hpp"
#include "details/pickle_additive.hpp"

namespace fwdpy11
{
    using MlocusAdditive
        = fwdpp_mlocus_gvalue<fwdpp::additive_diploid, pickle_additive>;

    inline MlocusAdditive
    create_MlocusAdditive(const fwdpp::trait t,
                          const fwdpy11::GeneticValueToFitnessMap& gv2w)
    {
        return MlocusAdditive(fwdpp::additive_diploid(t),
                              aggregate_additive_trait(), gv2w);
    }

    inline MlocusAdditive
    create_MlocusAdditive(const fwdpp::trait t,
                          const fwdpy11::GeneticValueToFitnessMap& gv2w,
                          const fwdpy11::GeneticValueNoise& noise)
    {
        return MlocusAdditive(fwdpp::additive_diploid(t),
                              aggregate_additive_trait(), gv2w, noise);
    }

    inline MlocusAdditive
    create_MlocusAdditive(const fwdpp::fitness t)
    {
        return MlocusAdditive(fwdpp::additive_diploid(t),
                              aggregate_additive_fitness());
    }

    inline MlocusAdditive
    create_MlocusAdditive(const fwdpp::fitness t,
                          const fwdpy11::GeneticValueToFitnessMap& gv2w,
                          const fwdpy11::GeneticValueNoise& noise)
    {
        return MlocusAdditive(fwdpp::additive_diploid(t),
                              aggregate_additive_fitness(), gv2w, noise);
    }
} // namespace fwdpy11
#endif
