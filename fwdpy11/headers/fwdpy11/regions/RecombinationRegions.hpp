#ifndef FWDPY11_RECOMBINATIONREGIONS_HPP
#define FWDPY11_RECOMBINATIONREGIONS_HPP

#include <limits>
#include <vector>
#include <algorithm>
#include <functional>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpy11/rng.hpp>
#include "Region.hpp"
#include "GeneticMapUnit.hpp"

namespace fwdpy11
{
    struct GeneticMap
    {
        virtual ~GeneticMap() = default;
        virtual std::vector<double> operator()(const GSLrng_t& rng) const = 0;
    };

    struct RecombinationRegions : public GeneticMap
    {
        std::vector<Region> regions;
        std::vector<double> weights;
        fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
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
                    lookup.reset(gsl_ran_discrete_preproc(weights.size(),
                                                          weights.data()));
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

    struct GeneralizedGeneticMap : public GeneticMap
    {
        std::vector<std::unique_ptr<GeneticMapUnit>> callbacks;
        GeneralizedGeneticMap(std::vector<std::unique_ptr<GeneticMapUnit>> c)
            : callbacks(std::move(c))
        {
        }

        std::vector<double>
        operator()(const GSLrng_t& rng) const final
        {
            std::vector<double> rv;
            for (auto&& c : callbacks)
                {
                    c->operator()(rng, rv);
                }
            if (!rv.empty())
                {
                    std::sort(begin(rv), end(rv));
                    rv.push_back(std::numeric_limits<double>::max());
                }
            return rv;
        }
    };

    struct MlocusRecombinationRegions
    {
        //TODO: rename
        std::vector<RecombinationRegions> regions;
        MlocusRecombinationRegions() : regions() {}
        inline void
        append(RecombinationRegions& r)
        {
            regions.emplace_back(std::move(r));
        }
    };
} // namespace fwdpy11

#endif
