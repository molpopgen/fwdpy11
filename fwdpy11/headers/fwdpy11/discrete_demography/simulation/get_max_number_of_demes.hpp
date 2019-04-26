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

#ifndef FWDPY11_GET_MAX_NUMBER_OF_DEMES_HPP
#define FWDPY11_GET_MAX_NUMBER_OF_DEMES_HPP

#include "../DiscreteDemography.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class get_max_number_of_demes
        {
          private:
            template <typename T>
            std::int32_t
            update_maxdeme_from_demography(std::int32_t m, const T& t)
            {
                for (auto&& i : t)
                    {
                        m = std::max(m, i.deme);
                    }
                return m;
            }

          public:
            template <typename METADATATYPE>
            std::int32_t
            operator()(const std::vector<METADATATYPE>& input_metadata,
                       const DiscreteDemography& demography)
            // TODO: has to update to deal with size of migration matrix.
            // TODO: check if deme IDs are contiguous for a simulation
            {
                std::int32_t maxdeme_from_metadata = -1;
                std::int32_t maxdeme_from_demography = -1;

                for (auto&& i : input_metadata)
                    {
                        if (i.deme < 0)
                            {
                                throw std::invalid_argument(
                                    "input deme labels must be non-negative");
                            }
                        maxdeme_from_metadata
                            = std::max(i.deme, maxdeme_from_metadata);
                    }
                for (auto&& i : demography.mass_migrations)
                    {
                        maxdeme_from_demography
                            = std::max(maxdeme_from_demography, i.source);
                        maxdeme_from_demography
                            = std::max(maxdeme_from_demography, i.destination);
                    }
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, demography.set_growth_rates);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, demography.set_deme_sizes);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, demography.set_selfing_rates);
                maxdeme_from_demography = update_maxdeme_from_demography(
                    maxdeme_from_demography, demography.set_migration_rates);
                auto temp
                    = std::max(maxdeme_from_metadata, maxdeme_from_demography)
                      + 1;
                if (demography.migmatrix == nullptr)
                    {
                        // There is no migration, so we are done
                        return temp;
                    }
                // NOTE: need to check >= 0 b/c we init everything to -1
                if (maxdeme_from_demography >= 0
                    && static_cast<std::size_t>(temp)
                           > demography.migmatrix->npops)
                    {
                        throw std::invalid_argument(
                            "MigrationMatrix contains too few demes");
                    }
                return std::max(temp, static_cast<std::int32_t>(
                                          demography.migmatrix->npops));
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
