#pragma once

#include <memory>

#include "GeneticValueIsTrait.hpp"
#include "GSSmo.hpp"
#include "MultivariateGSSmo.hpp"

namespace fwdpy11
{
    class GaussianStabilizingSelection : public GeneticValueIsTrait
    {
      private:
        std::shared_ptr<GeneticValueToFitnessMap> pimpl;

      public:
        explicit GaussianStabilizingSelection(const GSSmo &input)
            : GeneticValueIsTrait(input.total_dim), pimpl(input.clone())
        {
        }

        explicit GaussianStabilizingSelection(const MultivariateGSSmo &input)
            : GeneticValueIsTrait(input.total_dim), pimpl(input.clone())
        {
        }

        double
        operator()(const DiploidGeneticValueToFitnessData data) const final {
            return this->pimpl->operator()(data);
        }

        void
        update(const DiploidPopulation &pop) final
        {
            return this->pimpl->update(pop);
        }

        std::shared_ptr<GeneticValueToFitnessMap>
        clone() const final
        {
            return this->pimpl->clone();
        }
    };
}
