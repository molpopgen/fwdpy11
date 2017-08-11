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
#ifndef FWDPY11_MULTILOCUS_HPP__
#define FWDPY11_MULTILOCUS_HPP__

#include <functional>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct interlocus_rec
    {
        enum RECMODEL : int
        {
            BINOMIAL,
            POISSON
        };
        using mtype = std::underlying_type<RECMODEL>::type;
        const double param;
        RECMODEL m;
        interlocus_rec(const double param_, const mtype m_)
            : param(param_), m(static_cast<RECMODEL>(m_))
        {
        }

        std::function<unsigned(void)>
        callback(const fwdpy11::GSLrng_t& rng) const
        {
            if (m == RECMODEL::BINOMIAL)
                {
                    return std::bind(gsl_ran_binomial, rng.get(), param, 1);
                }
            return std::bind(gsl_ran_poisson, rng.get(), param);
        }

        mtype get_model() const
        {
            return static_cast<mtype>(m);
        }
    };
}
#endif
