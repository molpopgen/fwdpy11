#ifndef FWDPY11_MUTATIONREGIONS_HPP
#define FWDPY11_MUTATIONREGIONS_HPP

#include <vector>
#include <memory>
#include <fwdpp/gsl_discrete.hpp>
#include <gsl/gsl_randist.h>
#include "Region.hpp"
#include "Sregion.hpp"

namespace fwdpy11
{
    struct MutationRegions
    {
        std::vector<std::unique_ptr<Sregion>> regions;
        std::vector<double> weights;
        fwdpp::gsl_ran_discrete_t_ptr lookup;
        MutationRegions(std::vector<std::unique_ptr<Sregion>>&& r,
                        std::vector<double>&& w)
            : regions(std::move(r)), weights(std::move(w)),
              lookup(weights.empty() ? nullptr
                                     : gsl_ran_discrete_preproc(
                                         weights.size(), weights.data()))
        {
        }

        inline static MutationRegions
        create(double pneutral, std::vector<double>& nweights,
               std::vector<double>& sweights,
               std::vector<std::unique_ptr<Sregion>>& nregions,
               std::vector<std::unique_ptr<Sregion>>& sregions)
        {
            // Have to reweight input weights
            double sum_neutral
                = std::accumulate(begin(nweights), end(nweights), 0.0);
            double sum_selected
                = std::accumulate(begin(sweights), end(sweights), 0.0);

            std::transform(
                begin(nweights), end(nweights), begin(nweights),
                [sum_neutral](double d) { return d / sum_neutral; });
            std::transform(begin(nweights), end(nweights), begin(nweights),
                           [pneutral](double d) { return d * pneutral; });
            std::transform(
                begin(sweights), end(sweights), begin(sweights),
                [sum_selected](double d) { return d / sum_selected; });
            std::transform(
                begin(sweights), end(sweights), begin(sweights),
                [pneutral](double d) { return d * (1. - pneutral); });

            std::vector<double> weights(begin(nweights), end(nweights));
            weights.insert(end(weights), begin(sweights), end(sweights));

            std::vector<std::unique_ptr<Sregion>> combined(
                std::make_move_iterator(begin(nregions)),
                std::make_move_iterator(end(nregions)));
            // Refactored Jan 16, 2020 to use the 
            // range-based for instead of insert.
            // The previous didn't compile in
            // "extreme" debugging mode on GCC.
            for (auto& sr : sregions)
                {
                    combined.emplace_back(std::move(sr));
                }
            return MutationRegions(std::move(combined), std::move(weights));
        }
    };
} // namespace fwdpy11

#endif
