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
#include <fwdpp/gsl_discrete.hpp>
#include "../MigrationMatrix.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct migration_lookup
        {
            std::vector<fwdpp::gsl_ran_discrete_t_ptr> lookups;
            const bool null_migmatrix;
            migration_lookup(std::int32_t maxdemes, bool isnull)
                : lookups(maxdemes), null_migmatrix(isnull)
            {
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
