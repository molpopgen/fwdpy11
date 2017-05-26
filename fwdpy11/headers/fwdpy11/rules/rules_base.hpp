#ifndef FWDPY_RULES_BASE_HPP
#define FWDPY_RULES_BASE_HPP

#include "fwdpy11/fitness/fitness.hpp"
#include "fwdpy11/types.hpp"
#include <fwdpp/internal/gsl_discrete.hpp>
#include <stdexcept>
#include <vector>

namespace fwdpy11
{
    struct single_region_rules_base
    {
        std::vector<double> fitnesses;
        KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
        double wbar;
        single_region_rules_base()
            : fitnesses(std::vector<double>()),
              lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
              wbar(0.0)
        {
        }

        single_region_rules_base(single_region_rules_base &&) = default;

        single_region_rules_base(const single_region_rules_base &rhs)
            : fitnesses(rhs.fitnesses),
              lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
              wbar(rhs.wbar)
        {
            if (!fitnesses.empty())
                lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(fitnesses.size(), &fitnesses[0]));
        }

        virtual ~single_region_rules_base() {}

        virtual double w(singlepop_t &pop, const single_locus_fitness_fxn &ff)
            = 0;

        //! \brief Pick parent one
        virtual size_t
        pick1(const GSLrng_t &rng, const singlepop_t &) const
        {
            return gsl_ran_discrete(rng.get(), lookup.get());
        }

        //! \brief Pick parent 2.  Parent 1's data are passed along for models
        //! where that is relevant
        virtual size_t
        pick2(const GSLrng_t &rng, const singlepop_t &, const std::size_t p1,
              const double f) const
        {
            return (f == 1. || (f > 0. && gsl_rng_uniform(rng.get()) < f))
                       ? p1
                       : gsl_ran_discrete(rng.get(), lookup.get());
        }

        //! \brief Update some property of the offspring based on properties of
        //! the parents
        virtual void update(const GSLrng_t &rng, diploid_t &dip,
                            const singlepop_t &pop, const std::size_t p1,
                            const std::size_t p2)
            = 0;
    };
}

#endif
