#ifndef FWDPY11_RULES_QTRAIT_HPP__
#define FWDPY11_RULES_QTRAIT_HPP__

#include "fwdpy11/rules/rules_base.hpp"
#include <functional>
#include <cmath>
#include <gsl/gsl_sf_pow_int.h>
/*
  Custom "rules" policy for single-region simulations of traits
  under Gaussian stabilizing selection.

  This file is partly KT having fun, but also an effort to stress-test fwdpp's
  "experimental" API.

  The advantage of this struct is:
  1. An offspring has its G and E values automatically assigned.  This allows
  us to record the
  exact fitnesses/heritabilities used in the simulation over time.
*/

namespace fwdpy11
{
    namespace qtrait
    {
        struct qtrait_model_rules : public fwdpy11::single_region_rules_base
        {
            using base_t = fwdpy11::single_region_rules_base;
            std::function<double(double)> trait_to_fitness;
            std::function<double(const double g, const fwdpy11::diploid_t &,
                                 const fwdpy11::diploid_t &)>
                noise_function;
            double sigE;
            qtrait_model_rules(std::function<double(double)> t2f,
                               std::function<double(
                                   const double g, const fwdpy11::diploid_t &,
                                   const fwdpy11::diploid_t &)>
                                   noise) noexcept(false)
                : base_t(), trait_to_fitness(std::move(t2f)),
                  noise_function(std::move(noise))
            /*!
              Constructor throws std::runtime_error if params are not valid.
            */
            {
            }

            qtrait_model_rules(qtrait_model_rules &&) = default;

            qtrait_model_rules(const qtrait_model_rules &rhs)
                : base_t(rhs), sigE(rhs.sigE)
            {
            }

            virtual void
            w(const dipvector_t &diploids, gcont_t &gametes,
              const mcont_t &mutations)
            {
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

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            virtual void
            update(const gsl_rng *r, diploid_t &offspring,
                   const diploid_t &parent1, const diploid_t &parent2,
                   const gcont_t &gametes, const mcont_t &mutations,
                   const singlepop_fitness_fxn &ff) noexcept
            {
                offspring.g = ff(offspring, gametes, mutations);
                offspring.e = noise_function(offspring.g, parent1, parent2);
                offspring.w = trait_to_fitness(offspring.g + offspring.e);
                assert(std::isfinite(offspring.w));
                return;
            }
        };
    } // namespace qtrait
} // namespace fwdpy
#endif
