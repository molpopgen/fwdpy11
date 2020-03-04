//
// Copyright (C) 2029 Kevin Thornton <krthornt@uci.edu>
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
//
#ifndef FWDPY11_NO_NOISE_HPP
#define FWDPY11_NO_NOISE_HPP

#include "GeneticValueNoise.hpp"

namespace fwdpy11
{
    struct NoNoise : public GeneticValueNoise
    {
        double
        operator()(const GSLrng_t& /*rng*/,
                   const DiploidMetadata& /*offspring_metadata*/,
                   const std::size_t /*parent1*/,
                   const std::size_t /*parent2*/,
                   const DiploidPopulation& /*pop*/) const override
        {
            return 0.;
        }

        void
        update(const DiploidPopulation& /*pop*/) override
        {
        }

        std::unique_ptr<GeneticValueNoise>
        clone() const override
        {
            return std::unique_ptr<NoNoise>(new NoNoise());
        }

        pybind11::object
        pickle() const override
        {
            return pybind11::bytes("NoNoise");
        }

        static inline const NoNoise
        unpickle(pybind11::object& o)
        {
            auto s = o.cast<std::string>();
            if (s != "NoNoise")
                {
                    throw std::runtime_error("invalid object state");
                }
            return NoNoise();
        }
    };
} // namespace fwdpy11

#endif
