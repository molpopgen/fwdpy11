#ifndef FWDPY11_RECOMBINATIONREGIONS_HPP
#define FWDPY11_RECOMBINATIONREGIONS_HPP

#include <cstdlib>
#include <limits>
#include <vector>
#include <algorithm>
#include <functional>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpp/genetic_map/genetic_map_unit.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpy11/rng.hpp>
#include "Region.hpp"
#include "fwdpp/util/validators.hpp"
#include "fwdpy11/gsl/gsl_error_handler_wrapper.hpp"

namespace fwdpy11
{
    struct GeneticMap
    {
        virtual ~GeneticMap() = default;
        GeneticMap() = default;
        GeneticMap(const GeneticMap&) = delete;
        GeneticMap(GeneticMap&&) = default;
        GeneticMap& operator=(const GeneticMap&) = delete;
        GeneticMap& operator=(GeneticMap&&) = default;
        virtual std::vector<double> operator()(const GSLrng_t& rng) const = 0;
    };

    struct RecombinationRegions : public GeneticMap
    {
        std::vector<Region> regions;
        std::vector<double> weights;
        fwdpp::gsl_ran_discrete_t_ptr lookup;
        double recrate;
        RecombinationRegions(double rate, const std::vector<Region> r)
            : regions(std::move(r)), weights{}, lookup(nullptr), recrate(rate)
        {
            for (auto& reg : regions)
                {
                    weights.push_back(reg.weight);
                }
            if (!weights.empty())
                {
                    lookup.reset(
                        gsl_ran_discrete_preproc(weights.size(), weights.data()));
                }
            else if (rate > 0.0)
                {
                    throw std::invalid_argument(
                        "recombination rate > 0 incompatible with empty "
                        "regions");
                }
        }

        std::vector<double>
        operator()(const GSLrng_t& rng) const final
        {
            unsigned nbreaks = gsl_ran_poisson(rng.get(), recrate);
            if (nbreaks == 0)
                {
                    return {};
                }
            std::vector<double> rv;
            rv.reserve(nbreaks + 1);
            for (unsigned i = 0; i < nbreaks; ++i)
                {
                    std::size_t x = gsl_ran_discrete(rng.get(), lookup.get());
                    rv.push_back(regions[x](rng));
                }
            std::sort(begin(rv), end(rv));
            rv.push_back(std::numeric_limits<double>::max());
            return rv;
        }
    };

    struct PoissonCrossoverGenerator
    {
        // Must return exactly 1 value per call
        virtual void breakpoint(const GSLrng_t& rng, std::vector<double>& breakpoints)
            = 0;
        virtual double left() = 0;
        virtual double right() = 0;
        virtual double mean_number_xovers() = 0;
        virtual ~PoissonCrossoverGenerator() = default;
        virtual std::unique_ptr<PoissonCrossoverGenerator> ll_clone() = 0;
    };

    struct NonPoissonCrossoverGenerator
    {
        virtual void breakpoint(const GSLrng_t& rng, std::vector<double>& breakpoints)
            = 0;
        virtual double left() = 0;
        virtual double right() = 0;
        virtual ~NonPoissonCrossoverGenerator() = default;
        virtual std::unique_ptr<NonPoissonCrossoverGenerator> ll_clone() = 0;
    };

    struct GeneralizedGeneticMap : public GeneticMap
    {
        std::vector<std::unique_ptr<PoissonCrossoverGenerator>> poisson_callbacks;
        std::vector<std::unique_ptr<NonPoissonCrossoverGenerator>> non_poisson_callbacks;
        fwdpp::gsl_ran_discrete_t_ptr poisson_lookup;
        double sum_poisson_means;
        GeneralizedGeneticMap(
            std::vector<std::unique_ptr<PoissonCrossoverGenerator>> pc,
            std::vector<std::unique_ptr<NonPoissonCrossoverGenerator>> nc)
            : poisson_callbacks(std::move(pc)), non_poisson_callbacks(std::move(nc)),
              poisson_lookup(nullptr), sum_poisson_means(0.0)
        {
            std::vector<double> means;
            for (auto& i : poisson_callbacks)
                {
                    sum_poisson_means += i->mean_number_xovers();
                    means.push_back(i->mean_number_xovers());
                }
            if (!means.empty())
                {
                    poisson_lookup.reset(
                        gsl_ran_discrete_preproc(means.size(), means.data()));
                }
        }

        std::vector<double>
        operator()(const GSLrng_t& rng) const final
        {
            std::vector<double> rv;
            auto nc = gsl_ran_poisson(rng.get(), sum_poisson_means);
            for (unsigned i = 0; i < nc; ++i)
                {
                    auto region = gsl_ran_discrete(rng.get(), poisson_lookup.get());
                    poisson_callbacks[region]->breakpoint(rng, rv);
                }
            for (auto& i : non_poisson_callbacks)
                {
                    i->breakpoint(rng, rv);
                }
            if (!rv.empty())
                {
                    std::sort(begin(rv), end(rv));
                    rv.push_back(std::numeric_limits<double>::max());
                }
            return rv;
        }
    };
} // namespace fwdpy11

#endif
