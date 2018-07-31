#ifndef FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_MLOCUS_GVALUE_HPP__
#define FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_MLOCUS_GVALUE_HPP__

#include <vector>
#include <type_traits>
#include <functional>
#include <algorithm>
#include "../MlocusPopGeneticValueWithMapping.hpp"
#include "../noise.hpp"

namespace fwdpy11
{
    template <typename fwdppT, typename pickleFunction>
    struct fwdpp_mlocus_gvalue
        : public fwdpy11::MlocusPopGeneticValueWithMapping
    {
        static_assert(
            std::is_convertible<pickleFunction, std::function<pybind11::object(
                                                    const fwdppT&)>>::value,
            "pickling function must be convertible to "
            "std::function<pybind11::object(const fwdppT*)>");
        using gvalue_map_ptr
            = std::unique_ptr<fwdpy11::GeneticValueToFitnessMap>;
        const fwdppT gv;
        std::function<double(const std::vector<double>&)> agg;
        mutable std::vector<double> per_locus_genetic_values;
        const pickleFunction pickle_fxn;

        template <typename forwarded_fwdppT, typename agg_t>
        fwdpp_mlocus_gvalue(forwarded_fwdppT&& gv_, agg_t&& agg_)
            : MlocusPopGeneticValueWithMapping{ GeneticValueIsFitness() },
              gv{ std::forward<forwarded_fwdppT>(gv_) },
              agg{ std::forward<agg_t>(agg_) }, per_locus_genetic_values{},
              pickle_fxn{ pickleFunction() }
        {
        }

        template <typename forwarded_fwdppT, typename agg_t>
        fwdpp_mlocus_gvalue(forwarded_fwdppT&& gv_, agg_t&& agg_,
                            const GeneticValueToFitnessMap& gv2w_)
            : MlocusPopGeneticValueWithMapping{ gv2w_ },
              gv{ std::forward<forwarded_fwdppT>(gv_) },
              agg{ std::forward<agg_t>(agg_) }, per_locus_genetic_values{},
              pickle_fxn{ pickleFunction() }

        {
        }

        template <typename forwarded_fwdppT, typename agg_t>
        fwdpp_mlocus_gvalue(forwarded_fwdppT&& gv_, agg_t&& agg_,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : MlocusPopGeneticValueWithMapping{ gv2w_, noise_ },
              gv{ std::forward<forwarded_fwdppT>(gv_) },
              agg{ std::forward<agg_t>(agg_) }, per_locus_genetic_values{},
              pickle_fxn{ pickleFunction() }
        {
        }

        inline double
        operator()(const std::size_t diploid_index,
                   const fwdpy11::MlocusPop& pop) const
        {
            per_locus_genetic_values.clear();
            // Back inserter amoritizes memory allocations very quickly.
            std::transform(pop.diploids[diploid_index].begin(),
                           pop.diploids[diploid_index].end(),
                           std::back_inserter(per_locus_genetic_values),
                           [this, &pop](const DiploidGenotype& g) {
                               return gv(g, pop.gametes, pop.mutations);
                           });
            return agg(per_locus_genetic_values);
        }

        inline void
        update(const fwdpy11::MlocusPop& pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }
        virtual pybind11::object
        pickle() const
        {
            return pickle_fxn(gv);
        }
    };
} // namespace fwdpy11
#endif
