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

#ifndef FWDPY11_MIGRATION_LOOKUP_HPP
#define FWDPY11_MIGRATION_LOOKUP_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>
#include <fwdpp/gsl_discrete.hpp>
#include "../../rng.hpp"
#include "deme_property_types.hpp"
#include "../MigrationMatrix.hpp"
#include "../exceptions.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct migration_lookup
        {
            std::vector<fwdpp::gsl_ran_discrete_t_ptr> lookups, olookups;
            const bool null_migmatrix;
            migration_lookup(std::int32_t maxdemes, bool isnull)
                : lookups(maxdemes), olookups(maxdemes), null_migmatrix(isnull)
            {
            }
        };

        void
        build_migration_lookup(
            const std::unique_ptr<MigrationMatrix>& M,
            const current_deme_sizes_vector& current_deme_sizes,
            const selfing_rates_vector& selfing_rates, migration_lookup& ml)
        // The lookup data is the transpose of a migration matrix multiplied
        // by the source population sizes, element-wise per row.
        // NOTE: the migration matrix entered into a simulation could be transposed immediately, right?
        {
            if (M != nullptr)
                {
                    std::vector<double> temp;
                    std::size_t npops = ml.lookups.size();
                    temp.reserve(npops);
                    const auto& ref = current_deme_sizes.get();
                    for (std::size_t dest = 0; dest < npops; ++dest)
                        {
                            for (std::size_t source = 0; source < npops;
                                 ++source)
                                {
                                    // By default, input migration rates are
                                    // weighted by the current deme size...
                                    double scaling_factor
                                        = static_cast<double>(ref[source]);
                                    // ...unless we are told not to do that.
                                    // But if the deme size is zero, we make
                                    // sure it is removed as a possible source
                                    // of a parent.
                                    if (M->scaled == false && ref[source] == 0)
                                        {
                                            scaling_factor = 0.0;
                                        }
                                    temp.push_back(
                                        scaling_factor
                                        * M->M[source * npops + dest]);
                                }
                            if (std::find_if(begin(temp), end(temp),
                                             [](double d) { return d != 0.; })
                                != end(temp))
                                {
                                    ml.lookups[dest].reset(
                                        gsl_ran_discrete_preproc(temp.size(),
                                                                 temp.data()));
                                    bool nonzero = false;
                                    for (std::size_t source = 0;
                                         source < npops; ++source)
                                        {
                                            double p = 1.0
                                                       - selfing_rates
                                                             .get()[source];
                                            temp[source] *= p;
                                            if (p > 0.)
                                                {
                                                    nonzero = true;
                                                }
                                        }
                                    if (nonzero)
                                        {
                                            ml.olookups[dest].reset(
                                                gsl_ran_discrete_preproc(
                                                    temp.size(), temp.data()));
                                        }
                                    else
                                        {
                                            ml.olookups[dest].reset(nullptr);
                                        }
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
