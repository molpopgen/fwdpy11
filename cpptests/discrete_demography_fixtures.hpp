#pragma once

#include <memory>
#include <fwdpy11/rng.hpp>
#include "fwdpy11/discrete_demography/SetDemeSize.hpp"
#include "make_DiscreteDemography.hpp"
#include "discrete_demography_roundtrips.hpp"

struct population_fixture
{
    fwdpy11::GSLrng_t rng;
    fwdpy11::DiploidPopulation pop;
    std::vector<fwdpy11::discrete_demography::MassMigration> mass_migrations;
    std::vector<fwdpy11::discrete_demography::SetExponentialGrowth> set_growth_rates;
    std::vector<fwdpy11::discrete_demography::SetDemeSize> set_deme_sizes;
    std::vector<fwdpy11::discrete_demography::SetMigrationRates> set_migration_rates;
    std::vector<fwdpy11::discrete_demography::SetSelfingRate> set_selfing_rates;
    fwdpy11::discrete_demography::MigrationMatrix migmatrix;

    struct gsl_matrix_deleter
    {
        inline void
        operator()(gsl_matrix *m) const
        {
            gsl_matrix_free(m);
        }
    };

    using gsl_matrix_ptr = std::unique_ptr<gsl_matrix, gsl_matrix_deleter>;

    population_fixture()
        : rng{416134}, pop(100, 1.0), mass_migrations{}, set_growth_rates{},
          set_deme_sizes{}, set_migration_rates{}, set_selfing_rates{}, migmatrix{}
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
        migmatrix
            = fwdpy11::discrete_demography::MigrationMatrix(std::forward<Args>(args)...);
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
