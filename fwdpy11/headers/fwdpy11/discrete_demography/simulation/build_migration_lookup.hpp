//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
#ifndef FWDPY11_BUILD_MIGRATION_LOOKUP_HPP
#define FWDPY11_BUILD_MIGRATION_LOOKUP_HPP

#include <algorithm>
#include <sstream>
#include <functional>
#include "../../rng.hpp"
#include <memory>
#include "../MigrationMatrix.hpp"
#include "deme_property_types.hpp"
#include "../exceptions.hpp"

#include "migration_lookup.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        inline void
        build_migration_lookup(const MigrationMatrix& M,
                               const current_deme_sizes_vector& current_deme_sizes,
                               const next_deme_sizes_vector& next_deme_sizes,
                               migration_lookup& ml)
        {
            if (!M.empty())
                {
                    std::vector<double> temp;
                    std::size_t npops = ml.lookups.size();
                    temp.reserve(npops);
                    const auto& current_deme_sizes_ref = current_deme_sizes.get();
                    const auto& next_deme_sizes_ref = next_deme_sizes.get();
                    for (std::size_t dest = 0; dest < npops; ++dest)
                        {
                            for (std::size_t source = 0; source < npops; ++source)
                                {
                                    // By default, input migration rates are
                                    // weighted by the current deme size...
                                    double scaling_factor = static_cast<double>(
                                        current_deme_sizes_ref[source]);
                                    // ...unless we are told not to do that.
                                    // But if the deme size is zero, we make
                                    // sure it is removed as a possible source
                                    // of a parent.
                                    if (M.scaled == false
                                        && current_deme_sizes_ref[source] != 0)
                                        {
                                            scaling_factor = 1.0;
                                        }
                                    double rate_in = M.M[dest * npops + source];
                                    if (rate_in > 0.
                                        && (current_deme_sizes_ref[source] == 0
                                            || next_deme_sizes_ref[dest] == 0))
                                        {
                                            if (current_deme_sizes_ref[dest] != 0)
                                                {
                                                    std::ostringstream o;
                                                    o << "non-zero migration "
                                                         "rate from "
                                                         "empty parental deme "
                                                      << source
                                                      << " into destination deme "
                                                      << dest;
                                                    throw DemographyError(o.str());
                                                }
                                            if (next_deme_sizes_ref[source] != 0)
                                                {
                                                    std::ostringstream o;
                                                    o << "non-zero migration "
                                                         "rate into "
                                                         "empty destination "
                                                         "deme "
                                                      << dest << " from source deme "
                                                      << source;
                                                    throw DemographyError(o.str());
                                                }
                                            else // both are zero
                                                {
                                                    std::ostringstream o;
                                                    o << "non-zero migration "
                                                         "from empty parental "
                                                         "deme "
                                                      << source
                                                      << " into empty "
                                                         "destination deme "
                                                      << dest;

                                                    throw DemographyError(o.str());
                                                }
                                        }
                                    temp.push_back(scaling_factor * rate_in);
                                }
                            if (std::find_if(begin(temp), end(temp),
                                             [](double d) { return d != 0.; })
                                != end(temp))
                                {
                                    ml.lookups[dest].reset(gsl_ran_discrete_preproc(
                                        temp.size(), temp.data()));
                                }
                            else // There is no possible migration into this deme
                                {
                                    ml.lookups[dest].reset(nullptr);
                                }
                            temp.clear();
                        }
                }
        }
    } // namespace discrete_demography
} // namespace fwdpy11
#endif
