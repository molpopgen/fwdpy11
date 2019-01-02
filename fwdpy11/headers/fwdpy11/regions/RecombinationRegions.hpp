#ifndef FWDPY11_RECOMBINATIONREGIONS_HPP
#define FWDPY11_RECOMBINATIONREGIONS_HPP

#include <limits>
#include <vector>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpy11/rng.hpp>
#include "Region.hpp"

namespace fwdpy11
{
    struct RecombinationRegions
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
            lookup.reset(
                gsl_ran_discrete_preproc(weights.size(), weights.data()));
        }

        inline std::vector<double>
        operator()(const GSLrng_t& rng) const
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
            rv.push_back(std::numeric_limits<double>::max());
            return rv;
        }
    }; 
} // namespace fwdpy11

#endif
