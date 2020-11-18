#pragma once

#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/discrete_demography/simulation/demographic_model_state.hpp>

struct demographic_events
{
    fwdpy11::discrete_demography::DiscreteDemography::mass_migration_vector
        mass_migrations;
    fwdpy11::discrete_demography::DiscreteDemography::set_growth_rates_vector
        set_growth_rates;
    fwdpy11::discrete_demography::DiscreteDemography::set_deme_sizes_vector
        set_deme_sizes;
    fwdpy11::discrete_demography::DiscreteDemography::set_migration_rates_vector
        set_migration_rates;
    fwdpy11::discrete_demography::DiscreteDemography::set_selfing_rates_vector
        set_selfing_rates;
    std::unique_ptr<fwdpy11::discrete_demography::MigrationMatrix> migmatrix;

    demographic_events()
        : mass_migrations{}, set_growth_rates{}, set_deme_sizes{}, set_migration_rates{},
          set_selfing_rates{}, migmatrix{nullptr}
    {
    }
};

inline void
add_events(demographic_events&)
{
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           fwdpy11::discrete_demography::DiscreteDemography::mass_migration_vector
               mass_migrations,
           Args&&... args)
{
    if (!events.mass_migrations.empty())
        {
            throw std::invalid_argument("multiple mass migration event lists");
        }
    events.mass_migrations = std::move(mass_migrations);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           fwdpy11::discrete_demography::DiscreteDemography::set_deme_sizes_vector
               set_deme_sizes,
           Args&&... args)
{
    if (!events.set_deme_sizes.empty())
        {
            throw std::invalid_argument("multiple deme size change event lists");
        }
    events.set_deme_sizes = std::move(set_deme_sizes);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           fwdpy11::discrete_demography::DiscreteDemography::set_growth_rates_vector
               set_growth_rates,
           Args&&... args)
{
    if (!events.set_growth_rates.empty())
        {
            throw std::invalid_argument("multiple growth rate change event lists");
        }
    events.set_growth_rates = std::move(set_growth_rates);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           fwdpy11::discrete_demography::DiscreteDemography::set_migration_rates_vector
               set_migration_rates,
           Args&&... args)
{
    if (!events.set_migration_rates.empty())
        {
            throw std::invalid_argument("multiple migration rate change event lists");
        }
    events.set_migration_rates = std::move(set_migration_rates);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           fwdpy11::discrete_demography::DiscreteDemography::set_selfing_rates_vector
               set_selfing_rates,
           Args&&... args)
{
    if (!events.set_selfing_rates.empty())
        {
            throw std::invalid_argument("multiple selfing rate change event lists");
        }
    events.set_selfing_rates = std::move(set_selfing_rates);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline void
add_events(demographic_events& events,
           std::unique_ptr<fwdpy11::discrete_demography::MigrationMatrix> migmatrix,
           Args&&... args)
{
    if (events.migmatrix != nullptr)
        {
            throw std::invalid_argument("multiple migration matrices");
        }
    events.migmatrix = std::move(migmatrix);
    add_events(events, std::forward<Args>(args)...);
}

template <typename... Args>
inline fwdpy11::discrete_demography::DiscreteDemography
make_DiscreteDemography(Args&&... args)
{
    demographic_events events;
    add_events(events, std::forward<Args>(args)...);

    return fwdpy11::discrete_demography::DiscreteDemography(
        std::move(events.mass_migrations), std::move(events.set_growth_rates),
        std::move(events.set_deme_sizes), std::move(events.set_selfing_rates),
        std::move(events.migmatrix), std::move(events.set_migration_rates));
}
