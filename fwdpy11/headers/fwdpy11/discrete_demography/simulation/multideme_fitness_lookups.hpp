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

#ifndef FWDPY11_MULTIDEME_FITNESS_LOOKUPS_HPP
#define FWDPY11_MULTIDEME_FITNESS_LOOKUPS_HPP

#include <algorithm>
#include <numeric>
#include <limits>
#include "deme_property_types.hpp"
#include "../exceptions.hpp"
#include "../../rng.hpp"
#include <fwdpp/gsl_discrete.hpp>

namespace fwdpy11
{
    namespace discrete_demography
    {
        template <typename T> struct multideme_fitness_lookups
        {
            std::vector<T> starts, stops, offsets;
            std::vector<double> fitnesses;
            std::vector<std::uint32_t> individuals;
            std::vector<fwdpp::gsl_ran_discrete_t_ptr> lookups;

            multideme_fitness_lookups(std::int32_t max_number_demes)
                : starts(max_number_demes, 0), stops(max_number_demes, 0),
                  offsets(max_number_demes, 0), fitnesses(), individuals(),
                  lookups(max_number_demes)
            {
            }

            template <typename METADATATYPE>
            void
            update(const current_deme_sizes_vector& deme_sizes,
                   const std::vector<METADATATYPE>& metadata)
            {
                auto& deme_sizes_ref = deme_sizes.get();
                fitnesses.resize(std::accumulate(begin(deme_sizes_ref),
                                                 end(deme_sizes_ref), 0),
                                 -1.0);
                individuals.resize(fitnesses.size(),
                                   std::numeric_limits<std::uint32_t>::max());
                std::fill(begin(fitnesses), end(fitnesses),
                          -1); // TODO: remove
                std::partial_sum(begin(deme_sizes_ref), end(deme_sizes_ref),
                                 begin(stops));
                std::copy(begin(stops), end(stops) - 1, begin(starts) + 1);
                std::fill(begin(offsets), end(offsets), 0);
                for (auto&& md : metadata)
                    {
                        auto i = starts[md.deme] + offsets[md.deme];
                        fitnesses[i] = md.w;
                        individuals[i] = md.label;
                        offsets[md.deme]++;
                    }
                for (std::size_t i = 0; i < starts.size(); ++i)
                    {
                        if (deme_sizes_ref[i] > 0)
                            {
                                // NOTE: the size of the i-th deme's
                                // fitness array is starts[i]-stops[i]
                                lookups[i].reset(gsl_ran_discrete_preproc(
                                    stops[i] - starts[i],
                                    fitnesses.data() + starts[i]));
                            }
                    }
            }

            T
            get_parent(const GSLrng_t& rng,
                       const current_deme_sizes_vector& deme_sizes,
                       const std::int32_t deme) const
            {
                if (deme_sizes.get()[deme] == 0)
                    {
                        throw EmptyDeme("parental deme is empty");
                    }
                auto o = gsl_ran_discrete(rng.get(), lookups[deme].get());
                return individuals[starts[deme] + o];
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
