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
#ifndef FWDPY11_EVOLVETS_RECORDERS_HPP
#define FWDPY11_EVOLVETS_RECORDERS_HPP

#include <vector>
#include <algorithm>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <gsl/gsl_randist.h>
#include "SampleRecorder.hpp"

namespace fwdpy11
{
    struct no_ancient_samples
    /*! When no ancient samples are tracked, 
     * this will provide the most efficient 
     * way to "do nothing" b/c the 
     * repeated C++/Py round trips are avoided.
     */
    {
        inline void
        operator()(const DiploidPopulation&, SampleRecorder&) const
        {
        }
    };

    struct random_ancient_samples
    {
        GSLrng_t rng;
        std::vector<fwdpp::uint_t> individuals, timepoints;
        fwdpp::uint_t samplesize;
        std::size_t next_timepoint;

        random_ancient_samples(const std::uint32_t seed, const fwdpp::uint_t n,
                               std::vector<fwdpp::uint_t> t)
            : rng(seed), individuals{}, timepoints(std::move(t)),
              samplesize(n), next_timepoint(0)
        {
        }

        template <typename poptype>
        void
        sample(const poptype& pop, SampleRecorder& sr)
        {
            if (next_timepoint < timepoints.size()
                && pop.generation == timepoints[next_timepoint])
                {
                    individuals.resize(pop.N);
                    std::iota(individuals.begin(), individuals.end(), 0);
                    sr.samples.resize(samplesize);
                    gsl_ran_choose(rng.get(), sr.samples.data(), samplesize,
                                   individuals.data(), pop.N,
                                   sizeof(fwdpp::uint_t));
                    ++next_timepoint;
                }
        }

        inline void
        operator()(const DiploidPopulation& pop, SampleRecorder& sr)
        {
            sample(pop, sr);
        }
    };

} // namespace fwdpy11
#endif
