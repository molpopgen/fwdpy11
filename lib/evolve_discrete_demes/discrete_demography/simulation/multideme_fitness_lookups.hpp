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
#include <sstream>
#include <numeric>
#include <limits>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <fwdpy11/rng.hpp>
#include <core/gsl/gsl_discrete.hpp>
#include <fwdpp/gsl_discrete.hpp>
#include <stdexcept>
#include "multideme_fitness_bookmark.hpp"

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        template <typename T> struct multideme_fitness_lookups
        {
            std::vector<fwdpp::gsl_ran_discrete_t_ptr> lookups;

            multideme_fitness_lookups(std::int32_t max_number_demes)
                : lookups(max_number_demes)
            {
            }

            // NOTE: We require this to compile for some reason?
            //multideme_fitness_lookups(const multideme_fitness_lookups& other)
            //    : starts(other.starts), stops(other.stops), offsets(other.offsets),
            //      fitnesses(other.fitnesses), individuals(other.individuals),
            //      lookups(copy_lookups(other))
            //{
            //}

            //multideme_fitness_lookups&
            //operator=(const multideme_fitness_lookups& other)
            //{
            //    starts = other.starts;
            //    stops = other.stops;
            //    fitnesses = other.fitnesses;
            //    individuals = other.individuals;
            //    lookups = copy_lookups(other);
            //    return *this;
            //}

            //std::vector<fwdpp::gsl_ran_discrete_t_ptr>
            //copy_lookups(const multideme_fitness_lookups& other)
            //{
            //    std::vector<fwdpp::gsl_ran_discrete_t_ptr> rv;
            //    rv.resize(other.lookups.size());
            //    for (std::size_t i = 0; i < other.starts.size(); ++i)
            //        {
            //            if (other.stops[i] - other.starts[i] > 0)
            //                {
            //                    // NOTE: the size of the i-th deme's
            //                    // fitness array is starts[i]-stops[i]
            //                    rv[i].reset(gsl_ran_discrete_preproc(
            //                        other.stops[i] - other.starts[i],
            //                        other.fitnesses.data() + starts[i]));
            //                }
            //            else
            //                {
            //                    rv[i].reset(nullptr);
            //                }
            //        }

            //    return rv;
            //}

            void
            update(const multideme_fitness_bookmark& fitness_bookmark)
            {
                for (std::size_t i = 0; i < fitness_bookmark.starts.size(); ++i)
                    {
                        if (fitness_bookmark.stops[i] - fitness_bookmark.starts[i] > 0)
                            {
                                // NOTE: the size of the i-th deme's
                                // fitness array is starts[i]-stops[i]
                                //
                                // FIXME: this is the cause of the segfault
                                // in testing
                                auto first = fitness_bookmark.individual_fitness.data()
                                             + fitness_bookmark.starts[i];
                                auto last = first
                                            + (fitness_bookmark.stops[i]
                                               - fitness_bookmark.starts[i]

                                            );

                                // generate a more useful error message here...
                                if (std::all_of(first, last,
                                                [](const double d) { return d == 0.0; }))
                                    {
                                        std::ostringstream o;
                                        o << "all fitness values in deme " << i
                                          << " are zero.";
                                        throw fwdpy11::discrete_demography::
                                            LocalExtinction(o.str());
                                    }
                                // ...then use the call that skips the same check
                                fwdpy11_core::update_lookup_table_skip_zero_check(
                                    fitness_bookmark.individual_fitness.data()
                                        + fitness_bookmark.starts[i],
                                    fitness_bookmark.stops[i]
                                        - fitness_bookmark.starts[i],
                                    lookups[i]);
                            }
                        else
                            {
                                lookups[i].reset(nullptr);
                            }
                    }
            }

            T
            get_parent(const fwdpy11::GSLrng_t& rng,
                       const multideme_fitness_bookmark& fitness_bookmark,
                       const std::int32_t deme) const
            {
                auto o = gsl_ran_discrete(rng.get(), lookups[deme].get());
                return fitness_bookmark.individuals[fitness_bookmark.starts[deme] + o];
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11_core

#endif
