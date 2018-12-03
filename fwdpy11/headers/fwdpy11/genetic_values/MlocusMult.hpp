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
#ifndef FWDPY11_GENETIC_VALUES_MLOCUSMULT_HPP__
#define FWDPY11_GENETIC_VALUES_MLOCUSMULT_HPP__

#include "details/pickle_multiplicative.hpp"
#include "fwdpp_wrappers/fwdpp_mlocus_gvalue.hpp"

namespace fwdpy11
{
    using MlocusMult = fwdpp_mlocus_gvalue<fwdpp::multiplicative_diploid,
                                           pickle_multiplicative>;
    inline MlocusMult
    create_MlocusMult(const fwdpp::trait t,
                      const fwdpy11::GeneticValueToFitnessMap& gv2w)
    {
        return MlocusMult(fwdpp::multiplicative_diploid(t),
                          aggregate_mult_trait(), gv2w);
    }

    inline MlocusMult
    create_MlocusMult(const fwdpp::trait t,
                      const fwdpy11::GeneticValueToFitnessMap& gv2w,
                      const fwdpy11::GeneticValueNoise& noise)
    {
        return MlocusMult(fwdpp::multiplicative_diploid(t),
                          aggregate_mult_trait(), gv2w, noise);
    }

    inline MlocusMult
    create_MlocusMult(const fwdpp::fitness t)
    {
        return MlocusMult(fwdpp::multiplicative_diploid(t),
                          aggregate_mult_fitness());
    }

    inline MlocusMult
    create_MlocusMult(const fwdpp::fitness t,
                      const fwdpy11::GeneticValueToFitnessMap& gv2w,
                      const fwdpy11::GeneticValueNoise& noise)
    {
        return MlocusMult(fwdpp::multiplicative_diploid(t),
                          aggregate_mult_fitness(), gv2w, noise);
    }
} // namespace fwdpy11
#endif
