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
#ifndef FWDPY11_DEMOGRAPHIC_EVENTS_HPP
#define FWDPY11_DEMOGRAPHIC_EVENTS_HPP

#include <sstream>
#include <algorithm>
#include <memory>
#include <vector>
#include <gsl/gsl_matrix.h>

#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "simulation/detail.hpp"
#include "DiscreteDemographyState.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        class DiscreteDemography
        {
          private:
            template <typename T>
            void
            sort_events(std::vector<T>& v) noexcept
            {
                std::stable_sort(begin(v), end(v));
            }

            static std::string
            event_name(const SetDemeSize&)
            {
                return "SetDemeSize";
            }

            static std::string
            event_name(const SetExponentialGrowth&)
            {
                return "SetExponentialGrowth";
            }

            static std::string
            event_name(const SetSelfingRate&)
            {
                return "SetSelfingRate";
            }

            static std::string
            event_name(const SetMigrationRates&)
            {
                return "SetMigrationRates";
            }

            template <typename T>
            void
            validate_events(const std::vector<T>& v)
            {
                for (std::size_t i = 1; i < v.size(); ++i)
                    {
                        if (v[i].when == v[i - 1].when && v[i].deme == v[i - 1].deme)
                            {
                                std::ostringstream o;
                                o << "DiscreteDemography: multiple " << event_name(v[i])
                                  << " events the same deme in the same "
                                     "generation";
                                throw std::invalid_argument(o.str());
                            }
                    }
            }

            void
            validate_events(const std::vector<MassMigration>& v)
            // NOTE: this will have to be updated when we allow for
            // sex-specific migration.
            {
                for (std::size_t i = 1; i < v.size(); ++i)
                    {
                        if (v[i].when == v[i - 1].when && v[i].source == v[i - 1].source
                            && v[i].destination == v[i - 1].destination
                            && v[i].move_individuals == v[i - 1].move_individuals)
                            {
                                throw std::invalid_argument(
                                    "DiscreteDemography: multiple "
                                    "MassMigration "
                                    "events from the same source to the same "
                                    "destination in the same generation");
                            }
                    }
                // Check that all moves from deme i to elsewhere
                // at time t move <= 100% of deme i
                for (auto i = begin(v); i < end(v);)
                    {
                        if (i->move_individuals == true)
                            {
                                auto j = i + 1;
                                double sum = i->fraction;
                                for (; j < end(v) && j->source == i->source
                                       && i->when == j->when && j->move_individuals;
                                     ++j)
                                    {
                                        if (i->source == j->source)
                                            {
                                                sum += j->fraction;
                                            }
                                    }
                                if (sum > 1.0)
                                    {
                                        std::ostringstream o;
                                        o << "DiscreteDemography: at time " << i->when
                                          << ", attempting to move " << sum * 100.0
                                          << "% of deme " << i->source << " is invalid";
                                        throw std::invalid_argument(o.str());
                                    }
                                i = j;
                            }
                        else
                            {
                                ++i;
                            }
                    }
            }

            template <typename T>
            std::vector<T>
            init_events_vector(std::vector<T> v)
            {
                std::vector<T> rv(std::move(v));
                sort_events(rv);
                validate_events(rv);
                return rv;
            }

            template <typename T>
            void
            update_event_times(std::uint32_t t, T& range)
            {
                range.get().first = std::lower_bound(
                    range.get().first, range.get().second, t,
                    [](const typename T::value_type::first_type::value_type v,
                       std::uint32_t t) { return v.when < t; });
            }

            void
            check_if_no_migration()
            // If there are no nonzero off-diagonal elements,
            // and no migration rate changes during a sim,
            // then there is no migration. Thus, reset
            // migmatrix to empty.
            {
                if (migmatrix.empty())
                    {
                        return;
                    }
                if (!set_migration_rates.empty())
                    {
                        return;
                    }
                gsl_matrix_const_view v = gsl_matrix_const_view_array(
                    migmatrix.M.data(), migmatrix.npops, migmatrix.npops);
                gsl_vector_const_view diag = gsl_matrix_const_diagonal(&v.matrix);
                bool allequal = true;
                for (std::size_t i = 0; allequal == true && i < migmatrix.npops; ++i)
                    {
                        gsl_vector_const_view rv = gsl_matrix_const_row(&v.matrix, i);
                        double rsum = 0.0;
                        for (std::size_t j = 0; j < migmatrix.npops; ++j)
                            {
                                rsum += gsl_vector_get(&rv.vector, j);
                            }
                        if (rsum != gsl_vector_get(&diag.vector, i))
                            {
                                allequal = false;
                            }
                    }
                if (allequal)
                    {
                        migmatrix.M.clear();
                    }
            }

            void
            validate_change_migration_events()
            {
                if (migmatrix.empty())
                    {
                        if (!set_migration_rates.empty())
                            {
                                throw std::invalid_argument(
                                    "migration matrix is None but "
                                    "SetMigrationRates events are registered");
                            }
                        return;
                    }
                for (auto& event : set_migration_rates)
                    {
                        if (event.migrates.size() == migmatrix.npops)
                            {
                                if (migmatrix.scaled == true)
                                    {
                                        auto sum
                                            = std::accumulate(begin(event.migrates),
                                                              end(event.migrates), 0.);
                                        if (sum != 0.0 && sum != 1.)
                                            {
                                                throw std::invalid_argument(
                                                    "new migration rates must "
                                                    "sum to "
                                                    "1.0");
                                            }
                                    }
                            }
                        else if (event.migrates.size()
                                 != migmatrix.npops * migmatrix.npops)
                            {
                                throw std::invalid_argument("invalid matrix size");
                            }
                    }
            }

            DiscreteDemographyState model_state;
            std::vector<MassMigration> mass_migrations;
            std::vector<SetExponentialGrowth> set_growth_rates;
            std::vector<SetDemeSize> set_deme_sizes;
            std::vector<SetSelfingRate> set_selfing_rates;
            MigrationMatrix migmatrix;
            std::vector<SetMigrationRates> set_migration_rates;

          public:
            DiscreteDemography(std::vector<MassMigration> mass_migrations,
                               std::vector<SetExponentialGrowth> set_growth_rates,
                               std::vector<SetDemeSize> set_deme_sizes,
                               std::vector<SetSelfingRate> set_selfing_rates,
                               MigrationMatrix m,
                               std::vector<SetMigrationRates> set_migration_rates)
                : model_state(init_events_vector(mass_migrations),
                              init_events_vector(set_growth_rates),
                              init_events_vector(set_deme_sizes),
                              init_events_vector(set_selfing_rates), m,
                              init_events_vector(set_migration_rates)),
                  mass_migrations(init_events_vector(std::move(mass_migrations))),
                  set_growth_rates(init_events_vector(std::move(set_growth_rates))),
                  set_deme_sizes(init_events_vector(std::move(set_deme_sizes))),
                  set_selfing_rates(init_events_vector(std::move(set_selfing_rates))),
                  migmatrix(std::move(m)),
                  set_migration_rates(init_events_vector(std::move(set_migration_rates)))
            {
                check_if_no_migration();
                validate_change_migration_events();
            }

            void
            reset_model_state()
            {
                if (model_state.maxdemes > 0)
                    {
                        model_state = DiscreteDemographyState(
                            mass_migrations, set_growth_rates, set_deme_sizes,
                            set_selfing_rates, migmatrix, set_migration_rates);
                    }
            }

            const DiscreteDemographyState&
            get_model_state()
            // Not visible to Python
            {
                return model_state;
            }

            void
            set_model_state(DiscreteDemographyState state)
            // Not visible to Python
            {
                model_state = std::move(state);
            }

            void
            copy_state_to(DiscreteDemography& other)
            // Visible to Python as DiscreteDemography._clone_state_to
            {
                other.model_state = this->model_state;
            }

            // these getters are not exposed to Python

            auto
            get_mass_migrations() const
            {
                return mass_migrations;
            }

            auto
            get_set_deme_sizes() const
            {
                return set_deme_sizes;
            }

            auto
            get_set_selfing_rates() const
            {
                return set_selfing_rates;
            }

            auto
            get_set_growth_rates() const
            {
                return set_growth_rates;
            }

            auto
            get_set_migration_rates() const
            {
                return set_migration_rates;
            }

            auto
            get_migration_matrix() const
            {
                return migmatrix;
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
