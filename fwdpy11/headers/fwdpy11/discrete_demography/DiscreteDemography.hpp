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
#include <tuple>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <fwdpp/named_type.hpp>

#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "simulation/detail.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        using mass_migration_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<MassMigration>::const_iterator,
                      const std::vector<MassMigration>::const_iterator>,
            detail::mass_migration_range_t>;
        using deme_size_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetDemeSize>::const_iterator,
                      const std::vector<SetDemeSize>::const_iterator>,
            detail::deme_size_change_range_t>;
        using growth_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetExponentialGrowth>::const_iterator,
                      const std::vector<SetExponentialGrowth>::const_iterator>,
            detail::growth_rate_change_range_t>;
        using selfing_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetSelfingRate>::const_iterator,
                      const std::vector<SetSelfingRate>::const_iterator>,
            detail::selfing_rate_change_range_t>;

        using migration_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetMigrationRates>::const_iterator,
                      const std::vector<SetMigrationRates>::const_iterator>,
            detail::migration_rate_change_range_t>;

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
                        if (v[i].when == v[i - 1].when
                            && v[i].deme == v[i - 1].deme)
                            {
                                std::ostringstream o;
                                o << "DiscreteDemography: multiple "
                                  << event_name(v[i])
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
                        if (v[i].when == v[i - 1].when
                            && v[i].source == v[i - 1].source
                            && v[i].destination == v[i - 1].destination
                            && v[i].move_individuals
                                   == v[i - 1].move_individuals)
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
                                       && j->move_individuals;
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
                                        o << "DiscreteDemography: at time "
                                          << i->when << ", attempting to move "
                                          << sum << " of deme " << i->source
                                          << " is invalid";
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
            init_events_vector(std::vector<T>&& v)
            {
                std::vector<T> rv(std::move(v));
                sort_events(rv);
                validate_events(rv);
                return rv;
            }

            template <typename T>
            std::pair<typename std::vector<T>::const_iterator,
                      const typename std::vector<T>::const_iterator>
            set_range(const std::vector<T>& v)
            {
                return { v.cbegin(), v.cend() };
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
            // migmatrix to nullptr.
            {
                if (migmatrix == nullptr)
                    {
                        return;
                    }
                if (!set_migration_rates.empty())
                    {
                        return;
                    }
                gsl_matrix_const_view v = gsl_matrix_const_view_array(
                    migmatrix->M.data(), migmatrix->npops, migmatrix->npops);
                gsl_vector_const_view diag
                    = gsl_matrix_const_diagonal(&v.matrix);
                bool allequal = true;
                for (std::size_t i = 0;
                     allequal == true && i < migmatrix->npops; ++i)
                    {
                        gsl_vector_const_view rv
                            = gsl_matrix_const_row(&v.matrix, i);
                        double rsum = 0.0;
                        for (std::size_t j = 0; j < migmatrix->npops; ++j)
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
                        migmatrix.reset(nullptr);
                    }
            }

          public:
            using mass_migration_vector = std::vector<MassMigration>;
            using set_growth_rates_vector = std::vector<SetExponentialGrowth>;
            using set_deme_sizes_vector = std::vector<SetDemeSize>;
            using set_selfing_rates_vector = std::vector<SetSelfingRate>;
            using set_migration_rates_vector = std::vector<SetMigrationRates>;

            mass_migration_vector mass_migrations;
            set_growth_rates_vector set_growth_rates;
            set_deme_sizes_vector set_deme_sizes;
            set_selfing_rates_vector set_selfing_rates;
            std::unique_ptr<const MigrationMatrix> migmatrix;
            set_migration_rates_vector set_migration_rates;

            // pairs of iterators over the events
            mass_migration_range mass_migration_tracker;
            deme_size_change_range deme_size_change_tracker;
            growth_rate_change_range growth_rate_change_tracker;
            selfing_rate_change_range selfing_rate_change_tracker;
            migration_rate_change_range migration_rate_change_tracker;

            DiscreteDemography(mass_migration_vector mmig,
                               set_growth_rates_vector sg,
                               set_deme_sizes_vector size_changes,
                               set_selfing_rates_vector ssr,
                               std::unique_ptr<MigrationMatrix> m,
                               set_migration_rates_vector smr)
                : mass_migrations(init_events_vector(std::move(mmig))),
                  set_growth_rates(init_events_vector(std::move(sg))),
                  set_deme_sizes(init_events_vector(std::move(size_changes))),
                  set_selfing_rates(init_events_vector(std::move(ssr))),
                  migmatrix(std::move(m)),
                  set_migration_rates(init_events_vector(std::move(smr))),
                  mass_migration_tracker(set_range(mass_migrations)),
                  deme_size_change_tracker(set_range(set_deme_sizes)),
                  growth_rate_change_tracker(set_range(set_growth_rates)),
                  selfing_rate_change_tracker(set_range(set_selfing_rates)),
                  migration_rate_change_tracker(set_range(set_migration_rates))
            {
                check_if_no_migration();
            }

            void
            update_event_times(std::uint32_t current_pop_generation)
            // When a simulation starts with the population's generation time
            // not at zero, then we assume that the pop'n has been evolved
            // and we may need to update the iterators accordingly.
            // NOTE: needs test.
            {
                update_event_times(current_pop_generation,
                                   mass_migration_tracker);
                update_event_times(current_pop_generation,
                                   growth_rate_change_tracker);
                update_event_times(current_pop_generation,
                                   deme_size_change_tracker);
                update_event_times(current_pop_generation,
                                   selfing_rate_change_tracker);
                update_event_times(current_pop_generation,
                                   migration_rate_change_tracker);
            }
        };
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
