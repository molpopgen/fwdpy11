#ifndef FWDPY11_RULES_QTRAIT_HPP__
#define FWDPY11_RULES_QTRAIT_HPP__

#include "fwdpy11/rules/rules_base.hpp"
#include <fwdpy11/evolve/qtrait_api.hpp>
#include <pybind11/numpy.h>
#include <functional>
#include <cmath>
#include <gsl/gsl_sf_pow_int.h>

namespace fwdpy11
{
    namespace qtrait
    {
        struct qtrait_model_rules : public fwdpy11::single_region_rules_base
        {
            using base_t = fwdpy11::single_region_rules_base;
            trait_to_fitness_function trait_to_fitness;
            single_locus_noise_function noise_function;
            qtrait_model_rules(trait_to_fitness_function t2f,
                               single_locus_noise_function noise) noexcept
                : base_t(), trait_to_fitness(std::move(t2f)),
                  noise_function(std::move(noise))
            {
            }

            qtrait_model_rules(qtrait_model_rules &&) = default;

            qtrait_model_rules(const qtrait_model_rules &rhs) : base_t(rhs) {}

            virtual double
            w(singlepop_t &pop, const single_locus_fitness_fxn &ff)
            {
                auto N_curr = pop.diploids.size();
                if (fitnesses.size() < N_curr)
                    fitnesses.resize(N_curr);
                wbar = 0.;
                for (size_t i = 0; i < N_curr; ++i)
                    {
                        pop.diploids[i].g
                            = ff(pop.diploids[i], pop.gametes, pop.mutations);
                        pop.diploids[i].w = trait_to_fitness(
                            pop.diploids[i].g, pop.diploids[i].e);
                        assert(std::isfinite(pop.diploids[i].w));
                        fitnesses[i] = pop.diploids[i].w;
                        wbar += pop.diploids[i].w;
                    }
                wbar /= double(N_curr);
                lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
                return wbar;
            }

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            virtual void
            update(const GSLrng_t &rng, diploid_t &offspring,
                   const singlepop_t &pop, const std::size_t p1,
                   const std::size_t p2) noexcept
            {
                offspring.e
                    = noise_function(pop.diploids[p1], pop.diploids[p2]);
                return;
            }
        };

        struct qtrait_mloc_rules
        {
            mutable double wbar;

            multilocus_aggregator_function aggregator;
            trait_to_fitness_function trait_to_fitness;
            multilocus_noise_function noise_function;
            mutable std::vector<double> fitnesses;

            mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
            //! \brief Constructor
            qtrait_mloc_rules(multilocus_aggregator_function ag,
                              trait_to_fitness_function t2f,
                              multilocus_noise_function nf)
                : wbar(0.), aggregator{ std::move(ag) },
                  trait_to_fitness{ std::move(t2f) },
                  noise_function{ std::move(nf) }, fitnesses{}
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
            double
            w(multilocus_t &pop, const multilocus_genetic_value &gvalue) const
            {
                unsigned N_curr = pop.diploids.size();
                if (fitnesses.size() < N_curr)
                    fitnesses.resize(N_curr);
                wbar = 0.;

                for (unsigned i = 0; i < N_curr; ++i)
                    {
                        pop.diploids[i][0].g = aggregator(gvalue(
                            pop.diploids[i], pop.gametes, pop.mutations));
                        pop.diploids[i][0].w = trait_to_fitness(
                            pop.diploids[i][0].g, pop.diploids[i][0].e);
                        fitnesses[i] = pop.diploids[i][0].w;
                        wbar += fitnesses[i];
                    }

                wbar /= double(N_curr);

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
                return wbar;
            }

            //! \brief Pick parent one
            inline size_t
            pick1(const GSLrng_t &rng, const multilocus_t &pop) const
            {
                return gsl_ran_discrete(rng.get(), lookup.get());
            }

            //! \brief Pick parent 2.  Parent 1's data are passed along for
            //! models where that is relevant
            inline size_t
            pick2(const GSLrng_t &rng, const multilocus_t &,
                  const std::size_t p1, const double f) const
            {
                return ((f == 1.)
                        || (f > 0. && gsl_rng_uniform(rng.get()) < f))
                           ? p1
                           : gsl_ran_discrete(rng.get(), lookup.get());
            }

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            void
            update(const GSLrng_t &rng, multilocus_diploid_t &offspring,
                   const multilocus_t &pop, const std::size_t p1,
                   const std::size_t p2) const
            {
                offspring[0].e
                    = noise_function(pop.diploids[p1], pop.diploids[p2]);
            }
        };
    } // namespace qtrait
} // namespace fwdpy
#endif
