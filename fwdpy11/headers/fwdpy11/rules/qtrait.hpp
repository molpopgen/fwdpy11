#ifndef FWDPY11_RULES_QTRAIT_HPP__
#define FWDPY11_RULES_QTRAIT_HPP__

#include "fwdpy11/rules/rules_base.hpp"
#include <pybind11/numpy.h>
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

            qtrait_model_rules(const qtrait_model_rules &rhs) : base_t(rhs) {}

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
                   const single_locus_fitness_fxn &ff) noexcept
            {
                offspring.g = ff(offspring, gametes, mutations);
                offspring.e = noise_function(offspring.g, parent1, parent2);
                offspring.w = trait_to_fitness(offspring.g + offspring.e);
                assert(std::isfinite(offspring.w));
                return;
            }
        };

        struct qtrait_mloc_rules
        {
            mutable double wbar;
            using aggregator_signature
                = std::function<double(const pybind11::array_t<double>)>;
            using trait_to_fitness_signature = std::function<double(double)>;
            using noise_function_signature = std::function<double(
                const double g, const fwdpy11::multilocus_diploid_t &,
                const fwdpy11::multilocus_diploid_t &)>;

            aggregator_signature aggregator;
            trait_to_fitness_signature trait_to_fitness;
            noise_function_signature noise_function;
            mutable std::vector<double> fitnesses;

            mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
            //! \brief Constructor
            qtrait_mloc_rules(aggregator_signature ag,
                              trait_to_fitness_signature t2f,
                              noise_function_signature nf)
                : wbar(0.), aggregator{ std::move(ag) },
                  trait_to_fitness{ std::move(t2f) },
                  noise_function{ std::move(nf) },
				  fitnesses{}
            {
            }

            qtrait_mloc_rules(qtrait_mloc_rules &&) = default;

            qtrait_mloc_rules(const qtrait_mloc_rules &rhs)
            {
                if (!fitnesses.empty())
                    lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                        gsl_ran_discrete_preproc(fitnesses.size(),
                                                 &fitnesses[0]));
            }

            //! \brief The "fitness manager"
            void
            w(const multilocus_t::dipvector_t &diploids, gcont_t &gametes,
              const mcont_t &mutations) const
            {
                unsigned N_curr = diploids.size();
                if (fitnesses.size() < N_curr)
                    fitnesses.resize(N_curr);
                wbar = 0.;

                for (unsigned i = 0; i < N_curr; ++i)
                    {
                        for (auto region : diploids[i])
                            {
                                gametes[region.first].n
                                    = gametes[region.second].n = 0;
                            }

                        // the g/e/w fields will be populated via update()
                        fitnesses[i] = diploids[i][0].w;
                        wbar += fitnesses[i];
                    }

                wbar /= double(diploids.size());

                /*!
                  Black magic alert:
                  fwdpp_internal::gsl_ran_discrete_t_ptr contains a
                  std::unique_ptr wrapping the GSL pointer.
                  This type has its own deleter, which is convenient, because
                  operator= for unique_ptrs automagically calls the deleter
                  before assignment!
                  Details:
                  http://www.cplusplus.com/reference/memory/unique_ptr/operator=

                  This only works b/c the rhs of the expression below may be
                  treated as an rvalue reference.
                */
                lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
            }

            //! \brief Pick parent one
            inline size_t
            pick1(const gsl_rng *r) const
            {
                return gsl_ran_discrete(r, lookup.get());
            }

            //! \brief Pick parent 2.  Parent 1's data are passed along for
            //! models where that is relevant
            inline size_t
            pick2(const gsl_rng *r, const size_t p1, const double f,
                  const fwdpy11::multilocus_diploid_t &,
                  const fwdpy11::gcont_t &, const fwdpy11::mcont_t &) const
            {
                return ((f == 1.) || (f > 0. && gsl_rng_uniform(r) < f))
                           ? p1
                           : gsl_ran_discrete(r, lookup.get());
            }

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            void
            update(const gsl_rng *r, multilocus_diploid_t &offspring,
                   const multilocus_diploid_t &parent1,
                   const multilocus_diploid_t &parent2, const gcont_t &gametes,
                   const mcont_t &mutations,
                   const multilocus_genetic_value &gv) const
            {
                offspring[0].g = aggregator(gv(offspring, gametes, mutations));
                offspring[0].e
                    = noise_function(offspring[0].g, parent1, parent2);
                offspring[0].w
                    = trait_to_fitness(offspring[0].g + offspring[0].e);
            }
        };
    } // namespace qtrait
} // namespace fwdpy
#endif
