/*!
  \brief "Rules" class for the standard W-F model
*/
#ifndef FWDPY11_WF_RULES_HPP__
#define FWDPY11_WF_RULES_HPP__

#include <cassert>
#include "rules_base.hpp"
#include <fwdpp/fitness_models.hpp>

namespace fwdpy11
{
    struct wf_rules : public fwdpy11::single_region_rules_base
    /*!
      \note This type is very much under development.
    */
    {
        using base_t = fwdpy11::single_region_rules_base;
        wf_rules() noexcept(true) : base_t()
        /*!
          Constructor
        */
        {
        }

        wf_rules(wf_rules &&) = default;

        wf_rules(const wf_rules &rhs) : base_t(rhs)
        /*!
          Copy constructor
        */
        {
        }

        virtual double
        w(SlocusPop &pop, const single_locus_fitness_fxn &ff)
        {
            auto N_curr = pop.diploids.size();
            if (fitnesses.size() < N_curr)
                fitnesses.resize(N_curr);
            wbar = 0.;
            for (size_t i = 0; i < N_curr; ++i)
                {
                    pop.diploid_metadata[i].w = pop.diploid_metadata[i].g
                        = ff(pop.diploids[i], pop.gametes, pop.mutations);
                    assert(std::isfinite(pop.diploid_metadata[i].w));
                    fitnesses[i] = pop.diploid_metadata[i].w;
                    wbar += pop.diploid_metadata[i].w;
                }
            wbar /= double(N_curr);
            lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
                gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
            return wbar;
        }

        //! \brief Update some property of the offspring based on properties of
        //! the parents
        virtual void
        update(const GSLrng_t &rng, dip_metadata &offspring_metadata,
               const SlocusPop &pop, const std::size_t p1,
               const std::size_t p2) noexcept
        {
            offspring_metadata.e = 0.0;
            offspring_metadata.g = 0.0;
            offspring_metadata.parents[0]=p1;
            offspring_metadata.parents[1]=p2;
            return;
        }
    };
}

#endif
