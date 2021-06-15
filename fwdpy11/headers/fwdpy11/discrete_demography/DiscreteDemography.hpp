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
#include <fwdpp/util/named_type.hpp>

#include "MassMigration.hpp"
#include "MigrationMatrix.hpp"
#include "SetDemeSize.hpp"
#include "SetExponentialGrowth.hpp"
#include "SetSelfingRate.hpp"
#include "SetMigrationRates.hpp"
#include "simulation/detail.hpp"
#include "simulation/multideme_fitness_lookups.hpp"
#include "simulation/deme_properties.hpp"
#include "simulation/migration_lookup.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        using mass_migration_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<MassMigration>::const_iterator,
                      std::vector<MassMigration>::const_iterator>,
            detail::mass_migration_range_t>;
        using deme_size_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetDemeSize>::const_iterator,
                      std::vector<SetDemeSize>::const_iterator>,
            detail::deme_size_change_range_t>;
        using growth_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetExponentialGrowth>::const_iterator,
                      std::vector<SetExponentialGrowth>::const_iterator>,
            detail::growth_rate_change_range_t>;
        using selfing_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetSelfingRate>::const_iterator,
                      std::vector<SetSelfingRate>::const_iterator>,
            detail::selfing_rate_change_range_t>;

        using migration_rate_change_range = fwdpp::strong_types::named_type<
            std::pair<std::vector<SetMigrationRates>::const_iterator,
                      std::vector<SetMigrationRates>::const_iterator>,
            detail::migration_rate_change_range_t>;

        class demographic_model_state
        /// Added in 0.6.0 to hold and manage
        /// the relevant data structures.
        {
          private:
            std::unique_ptr<MigrationMatrix>
            init_migmatrix(const std::unique_ptr<const MigrationMatrix>& Minput)
            {
                if (Minput == nullptr)
                    {
                        return nullptr;
                    }
                return std::unique_ptr<MigrationMatrix>(new MigrationMatrix(*Minput));
            }

            template <typename METADATATYPE>
            std::int32_t
            update_max_demes(const std::vector<METADATATYPE>& metadata,
                             const std::unique_ptr<const MigrationMatrix>& migmatrix,
                             std::int32_t maxdemes_from_demography)
            {
                std::int32_t max_from_metadata = -1;
                for (auto& md : metadata)
                    {
                        max_from_metadata = std::max(md.deme, max_from_metadata);
                    }
                auto temp = std::max(max_from_metadata + 1, maxdemes_from_demography);
                if (migmatrix == nullptr)
                    {
                        return temp;
                    }
                if (static_cast<std::size_t>(temp) > migmatrix->npops)
                    {
                        throw std::invalid_argument(
                            "MigrationMatrix contains too few demes");
                    }
                if (static_cast<std::size_t>(temp) < migmatrix->npops)
                    {
                        throw std::invalid_argument(
                            "MigrationMatrix contains too many demes");
                    }
                return std::max(temp, static_cast<std::int32_t>(migmatrix->npops));
            }

            std::uint32_t next_global_N;

          public:
            std::int32_t maxdemes;
            multideme_fitness_lookups<std::uint32_t> fitnesses;
            deme_properties sizes_rates;
            std::unique_ptr<MigrationMatrix> M;
            migration_lookup miglookup;

            // NOTE: demography.update_event_times() needs to have been
            // called first!
            template <typename METADATATYPE>
            demographic_model_state(
                const std::vector<METADATATYPE>& metadata,
                const std::unique_ptr<const MigrationMatrix>& migmatrix,
                std::int32_t maxdemes_from_demography)
                : next_global_N(0), maxdemes(update_max_demes(metadata, migmatrix,
                                                              maxdemes_from_demography)),
                  fitnesses(maxdemes), sizes_rates(maxdemes, metadata),
                  M(init_migmatrix(migmatrix)), miglookup(maxdemes, M == nullptr)
            {
            }

            // This constructor is only used when resetting
            // the state from an event like pickling a DiscreteDemography
            // instance.
            demographic_model_state(std::int32_t maxdemes_, deme_properties sizes_rates_,
                                    std::unique_ptr<MigrationMatrix> M_)
                : next_global_N(0), maxdemes(maxdemes_), fitnesses(maxdemes),
                  sizes_rates(std::move(sizes_rates_)), M(std::move(M_)),
                  miglookup(maxdemes, M == nullptr)
            {
                next_global_N
                    = std::accumulate(begin(sizes_rates.next_deme_sizes.get()),
                                      end(sizes_rates.next_deme_sizes.get()), 0u);
            }

            demographic_model_state(const demographic_model_state& other)
                : next_global_N{other.next_global_N}, maxdemes{other.maxdemes},
                  fitnesses{maxdemes}, sizes_rates{other.sizes_rates},
                  M{other.M == nullptr ? nullptr : new MigrationMatrix(*other.M)},
                  miglookup(maxdemes, M == nullptr)
            {
            }

            demographic_model_state(demographic_model_state&& other)
                : next_global_N{other.next_global_N}, maxdemes{other.maxdemes},
                  fitnesses{std::move(other.fitnesses)},
                  sizes_rates{std::move(other.sizes_rates)}, M{std::exchange(other.M,
                                                                             nullptr)},
                  miglookup(maxdemes, M == nullptr)
            {
            }

            demographic_model_state&
            operator=(const demographic_model_state& other)
            {
                next_global_N = other.next_global_N;
                maxdemes = other.maxdemes;
                fitnesses = multideme_fitness_lookups<std::uint32_t>{maxdemes};
                if (other.M == nullptr)
                    {
                        M = nullptr;
                    }
                else
                    {
                        M = std::make_unique<MigrationMatrix>(*other.M);
                    }
                auto lookup = migration_lookup(maxdemes, M == nullptr);
                miglookup.lookups.swap(lookup.lookups);
                miglookup.null_migmatrix = lookup.null_migmatrix;
                return *this;
            }

            demographic_model_state&
            operator=(demographic_model_state&& other)
            {
                next_global_N = other.next_global_N;
                maxdemes = other.maxdemes;
                fitnesses = std::move(other.fitnesses);
                M.swap(other.M);
                std::swap(miglookup, other.miglookup);
                return *this;
            }

            void
            set_next_global_N(std::uint32_t N)
            {
                next_global_N = N;
            }

            std::uint32_t
            ttlN_next() const
            {
                return next_global_N;
            }

            bool
            will_go_globally_extinct() const
            {
                return ttlN_next() == 0;
            }
        };

        using demographic_model_state_pointer = std::unique_ptr<demographic_model_state>;

        class DiscreteDemography
        {
          private:
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
                template <
                    typename mass_migration_vector, typename set_growth_rates_vector,
                    typename set_deme_sizes_vector, typename set_selfing_rates_vector,
                    typename set_migration_rates_vector>
                std::int32_t
                operator()(const std::unique_ptr<const MigrationMatrix>& migmatrix,
                           const mass_migration_vector& mass_migrations,
                           const set_growth_rates_vector& set_growth_rates,
                           const set_deme_sizes_vector& set_deme_sizes,
                           const set_selfing_rates_vector& set_selfing_rates,
                           const set_migration_rates_vector& set_migration_rates)
                // TODO: has to update to deal with size of migration matrix.
                // TODO: check if deme IDs are contiguous for a simulation
                {
                    std::int32_t maxdeme_from_demography = -1;

                    for (auto&& i : mass_migrations)
                        {
                            maxdeme_from_demography
                                = std::max(maxdeme_from_demography, i.source);
                            maxdeme_from_demography
                                = std::max(maxdeme_from_demography, i.destination);
                        }
                    maxdeme_from_demography = update_maxdeme_from_demography(
                        maxdeme_from_demography, set_growth_rates);
                    maxdeme_from_demography = update_maxdeme_from_demography(
                        maxdeme_from_demography, set_deme_sizes);
                    maxdeme_from_demography = update_maxdeme_from_demography(
                        maxdeme_from_demography, set_selfing_rates);
                    maxdeme_from_demography = update_maxdeme_from_demography(
                        maxdeme_from_demography, set_migration_rates);
                    auto temp = maxdeme_from_demography + 1;
                    if (migmatrix == nullptr)
                        {
                            // There is no migration, so we are done
                            return temp;
                        }
                    if (static_cast<std::size_t>(temp) > migmatrix->npops)
                        {
                            throw std::invalid_argument(
                                "MigrationMatrix contains too few demes");
                        }
                    return temp;
                }
            };

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
            init_events_vector(std::vector<T>&& v)
            {
                std::vector<T> rv(std::move(v));
                sort_events(rv);
                validate_events(rv);
                return rv;
            }

            template <typename T>
            std::pair<typename std::vector<T>::const_iterator,
                      typename std::vector<T>::const_iterator>
            set_range(const std::vector<T>& v)
            {
                return {v.cbegin(), v.cend()};
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

            template <typename T1, typename T2>
            void
            reset_range(T1& range, T2& events)
            {
                range.get().first = events.cbegin();
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
                gsl_vector_const_view diag = gsl_matrix_const_diagonal(&v.matrix);
                bool allequal = true;
                for (std::size_t i = 0; allequal == true && i < migmatrix->npops; ++i)
                    {
                        gsl_vector_const_view rv = gsl_matrix_const_row(&v.matrix, i);
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

            void
            validate_change_migration_events()
            {
                if (migmatrix == nullptr)
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
                        if (event.migrates.size() == migmatrix->npops)
                            {
                                if (migmatrix->scaled == true)
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
                                 != migmatrix->npops * migmatrix->npops)
                            {
                                throw std::invalid_argument("invalid matrix size");
                            }
                    }
            }

            demographic_model_state_pointer model_state;

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

            std::int32_t maxdemes_from_demographic_events;

            DiscreteDemography(mass_migration_vector mmig, set_growth_rates_vector sg,
                               set_deme_sizes_vector size_changes,
                               set_selfing_rates_vector ssr,
                               std::unique_ptr<MigrationMatrix> m,
                               set_migration_rates_vector smr)
                : model_state(nullptr),
                  mass_migrations(init_events_vector(std::move(mmig))),
                  set_growth_rates(init_events_vector(std::move(sg))),
                  set_deme_sizes(init_events_vector(std::move(size_changes))),
                  set_selfing_rates(init_events_vector(std::move(ssr))),
                  migmatrix(std::move(m)),
                  set_migration_rates(init_events_vector(std::move(smr))),
                  mass_migration_tracker(set_range(mass_migrations)),
                  deme_size_change_tracker(set_range(set_deme_sizes)),
                  growth_rate_change_tracker(set_range(set_growth_rates)),
                  selfing_rate_change_tracker(set_range(set_selfing_rates)),
                  migration_rate_change_tracker(set_range(set_migration_rates)),
                  maxdemes_from_demographic_events{get_max_number_of_demes()(
                      migmatrix, mass_migrations, set_growth_rates, set_deme_sizes,
                      set_selfing_rates, set_migration_rates)}
            {
                check_if_no_migration();
                validate_change_migration_events();
            }

            DiscreteDemography(const DiscreteDemography& other)
                : model_state{other.model_state == nullptr
                                  ? nullptr
                                  : new demographic_model_state{*other.model_state}},
                  mass_migrations(other.mass_migrations),
                  set_growth_rates(other.set_growth_rates),
                  set_deme_sizes(other.set_deme_sizes),
                  set_selfing_rates(other.set_selfing_rates),
                  migmatrix{other.migmatrix == nullptr
                                ? nullptr
                                : new MigrationMatrix(*other.migmatrix)},
                  set_migration_rates(other.set_migration_rates),
                  mass_migration_tracker(set_range(mass_migrations)),
                  deme_size_change_tracker(set_range(set_deme_sizes)),
                  growth_rate_change_tracker(set_range(set_growth_rates)),
                  selfing_rate_change_tracker(set_range(set_selfing_rates)),
                  migration_rate_change_tracker(set_range(set_migration_rates)),
                  maxdemes_from_demographic_events{
                      other.maxdemes_from_demographic_events}
            {
            }

            DiscreteDemography(DiscreteDemography&& other)
                : model_state{std::exchange(other.model_state, nullptr)},
                  mass_migrations(std::exchange(other.mass_migrations, {})),
                  set_growth_rates(std::exchange(other.set_growth_rates, {})),
                  set_deme_sizes(std::exchange(other.set_deme_sizes, {})),
                  set_selfing_rates(other.set_selfing_rates), migmatrix{std::exchange(
                                                                  other.migmatrix,
                                                                  nullptr)},
                  set_migration_rates(std::exchange(other.set_migration_rates, {})),
                  mass_migration_tracker(other.mass_migration_tracker),
                  deme_size_change_tracker(other.deme_size_change_tracker),
                  growth_rate_change_tracker(other.growth_rate_change_tracker),
                  selfing_rate_change_tracker(other.selfing_rate_change_tracker),
                  migration_rate_change_tracker(other.migration_rate_change_tracker),
                  maxdemes_from_demographic_events{
                      other.maxdemes_from_demographic_events}
            {
            }

            DiscreteDemography&
            operator=(const DiscreteDemography& other)
            {
                if (other.migmatrix == nullptr)
                    {
                        migmatrix = nullptr;
                    }
                else
                    {
                        migmatrix = std::make_unique<MigrationMatrix>(*other.migmatrix);
                    }
                mass_migrations = other.mass_migrations;
                set_growth_rates = other.set_growth_rates;
                set_deme_sizes = other.set_deme_sizes;
                set_selfing_rates = other.set_selfing_rates;

                mass_migration_tracker.get() = set_range(mass_migrations);
                deme_size_change_tracker.get() = set_range(set_deme_sizes);
                growth_rate_change_tracker.get() = set_range(set_growth_rates);
                selfing_rate_change_tracker.get() = set_range(set_selfing_rates);
                migration_rate_change_tracker.get() = set_range(set_migration_rates);

                maxdemes_from_demographic_events
                    = other.maxdemes_from_demographic_events;

                return *this;
            }

            DiscreteDemography&
            operator=(DiscreteDemography&& other)
            {
                migmatrix = std::exchange(other.migmatrix, nullptr);
                mass_migrations.swap(other.mass_migrations);
                set_growth_rates.swap(other.set_growth_rates);
                set_deme_sizes.swap(other.set_deme_sizes);
                set_selfing_rates.swap(other.set_selfing_rates);

                mass_migration_tracker.get() = set_range(mass_migrations);
                deme_size_change_tracker.get() = set_range(set_deme_sizes);
                growth_rate_change_tracker.get() = set_range(set_growth_rates);
                selfing_rate_change_tracker.get() = set_range(set_selfing_rates);
                migration_rate_change_tracker.get() = set_range(set_migration_rates);

                maxdemes_from_demographic_events
                    = other.maxdemes_from_demographic_events;

                return *this;
            }

            void
            update_event_times(std::uint32_t current_pop_generation)
            // When a simulation starts with the population's generation time
            // not at zero, then we assume that the pop'n has been evolved
            // and we may need to update the iterators accordingly.
            // NOTE: needs test.
            {
                reset_range(mass_migration_tracker, mass_migrations);
                update_event_times(current_pop_generation, mass_migration_tracker);
                reset_range(growth_rate_change_tracker, set_growth_rates);
                update_event_times(current_pop_generation, growth_rate_change_tracker);
                reset_range(deme_size_change_tracker, set_deme_sizes);
                update_event_times(current_pop_generation, deme_size_change_tracker);
                reset_range(selfing_rate_change_tracker, set_selfing_rates);
                update_event_times(current_pop_generation, selfing_rate_change_tracker);
                reset_range(migration_rate_change_tracker, set_migration_rates);
                update_event_times(current_pop_generation,
                                   migration_rate_change_tracker);
            }

            demographic_model_state_pointer
            get_model_state()
            // Not visible to Python
            {
                demographic_model_state_pointer rv(std::move(model_state));
                return rv;
            }

            void
            set_model_state(demographic_model_state_pointer state)
            // Not visible to Python
            {
                model_state = std::move(state);
            }
        };

        template <typename METADATATYPE>
        inline demographic_model_state_pointer
        initialize_model_state(std::uint32_t generation,
                               const std::vector<METADATATYPE>& metadata,
                               DiscreteDemography& demography)
        {
            // "steal" pointer from input
            auto rv = demography.get_model_state();

            if (rv == nullptr || generation == 0)
                // If there is no state, then we need to make
                // one.  If there is a state, but the generation
                // is zero, then we assume that the demography
                // has been used for a different simulatin replicate
                // and thus reset it.
                {
                    demography.update_event_times(generation);
                    rv.reset(new demographic_model_state(
                        metadata, demography.migmatrix,
                        demography.maxdemes_from_demographic_events));
                }
            return rv;
        }

        inline void
        save_model_state(demographic_model_state_pointer state,
                         DiscreteDemography& demography)
        {
            demography.set_model_state(std::move(state));
        }
    } // namespace discrete_demography
} // namespace fwdpy11

#endif
