#ifndef FWDPY11_RECOMBINATIONREGIONS_HPP
#define FWDPY11_RECOMBINATIONREGIONS_HPP

#include <vector>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <gsl/gsl_randist.h>
#include "Region.hpp"

namespace fwdpy11
{
    struct RecombinationRegions
    {
        std::vector<Region> regions;
        std::vector<double> weights;
        fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
        RecombinationRegions(const std::vector<Region> r)
            : regions(std::move(r)), weights{}, lookup(nullptr)
        {
            for (auto& reg : regions)
                {
                    weights.push_back(reg.weight);
                }
            lookup.reset(
                gsl_ran_discrete_preproc(weights.size(), weights.data()));
        }
    };
} // namespace fwdpy11

#endif
