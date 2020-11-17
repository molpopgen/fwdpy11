#pragma once

#include <memory>
#include <fwdpy11/rng.hpp>
#include "make_DiscreteDemography.hpp"
#include "discrete_demography_roundtrips.hpp"

struct mock_population_fixture
{
    fwdpy11::GSLrng_t rng;
    MockPopulation pop;
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

    struct gsl_matrix_deleter
    {
        inline void
        operator()(gsl_matrix *m) const
        {
            gsl_matrix_free(m);
        }
    };

    using gsl_matrix_ptr = std::unique_ptr<gsl_matrix, gsl_matrix_deleter>;

    mock_population_fixture()
        : rng{416134}, pop(100), mass_migrations{}, set_growth_rates{}, set_deme_sizes{},
          set_migration_rates{}, set_selfing_rates{}, migmatrix{nullptr}
    {
    }

    inline auto
    make_model()
    {
        return make_DiscreteDemography(
            std::move(mass_migrations), std::move(set_growth_rates),
            std::move(set_deme_sizes), std::move(set_selfing_rates),
            std::move(set_migration_rates), std::move(migmatrix));
    }

    template <typename... Args>
    inline void
    set_migmatrix(Args &&...args)
    {
        migmatrix = std::make_unique<fwdpy11::discrete_demography::MigrationMatrix>(
            std::forward<Args>(args)...);
    }

    auto
    make_gsl_matrix(std::size_t n)
    {
        gsl_matrix_ptr rv(gsl_matrix_calloc(n, n), gsl_matrix_deleter());
        return rv;
    }

    std::vector<double>
    convert_matrix(const gsl_matrix_ptr &mptr)
    {
        return std::vector<double>(mptr->data, mptr->data + mptr->size1 * mptr->size2);
    }
};
