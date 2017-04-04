/*!
  \brief "Rules" class for the standard W-F model
*/
#ifndef FWDPY11_WF_RULES_HPP__
#define FWDPY11_WF_RULES_HPP__

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

        virtual void
        w(const dipvector_t &diploids, gcont_t &gametes,
          const mcont_t &mutations)
        {
            index = 0; // reset this variable
            auto N_curr = diploids.size();
            if (fitnesses.size() < N_curr)
                fitnesses.resize(N_curr);
            wbar = 0.;
            for (size_t i = 0; i < N_curr; ++i)
                {
                    gametes[diploids[i].first].n
                        = gametes[diploids[i].second].n = 0;
                    fitnesses[i] = diploids[i].w;
                    wbar += diploids[i].w;
                }
            wbar /= double(N_curr);
            lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
        }

        //! \brief Update some property of the offspring based on properties of
        //! the parents
        virtual void
        update(const gsl_rng *r, diploid_t &offspring, const diploid_t &,
               const diploid_t &, const gcont_t &gametes,
               const mcont_t &mutations,
               const singlepop_fitness_fxn &ff) noexcept
        {
            offspring.w = ff(offspring, gametes, mutations);
            offspring.e = 0.0;
            offspring.g = 0.0;
            offspring.label = index++;
            assert(std::isfinite(offspring.w));
            return;
        }
    };
}

#endif
