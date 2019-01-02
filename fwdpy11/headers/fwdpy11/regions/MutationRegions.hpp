#ifndef FWDPY11_MUTATIONREGIONS_HPP
#define FWDPY11_MUTATIONREGIONS_HPP

#include <vector>
#include <memory>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <gsl/gsl_randist.h>
#include "Region.hpp"
#include "Sregion.hpp"

namespace fwdpy11
{
    struct MutationRegions
    {
        std::vector<std::unique_ptr<Sregion>> regions;
        std::vector<double> weights;
        fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
        MutationRegions(std::vector<std::unique_ptr<Sregion>>&& r,
                        std::vector<double>&& w)
            : regions(std::move(r)), weights(std::move(w)),
              lookup(gsl_ran_discrete_preproc(weights.size(), weights.data()))
        {
        }
    };
} // namespace fwdpy11

#endif
