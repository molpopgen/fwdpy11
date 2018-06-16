#ifndef FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_SLOCUS_GVALUE_HPP__
#define FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_SLOCUS_GVALUE_HPP__

#include "../SlocusPopGeneticValueWithMapping.hpp"
#include "../noise.hpp"

namespace fwdpy11
{
    template <typename fwdppT>
    struct fwdpp_slocus_gvalue
        : public fwdpy11::SlocusPopGeneticValueWithMapping
    {
        using gvalue_map_ptr
            = std::unique_ptr<fwdpy11::GeneticValueToFitnessMap>;
        const fwdppT gv;

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_)
            : SlocusPopGeneticValueWithMapping{ std::unique_ptr<
                  GeneticValueToFitnessMap>(new GeneticValueIsFitness()) },
              gv{ std::forward<forwarded_fwdppT>(gv_) }
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_)
            : SlocusPopGeneticValueWithMapping{ gv2w_.clone() }, gv{
                  std::forward<forwarded_fwdppT>(gv_)
              }
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : SlocusPopGeneticValueWithMapping{ gv2w_.clone(),
                                                noise_.clone() },
              gv{ std::forward<forwarded_fwdppT>(gv_) }

        {
        }

        inline double
        operator()(const std::size_t diploid_index,
                   const fwdpy11::SlocusPop& pop) const
        {
            return gv(pop.diploids[diploid_index], pop.gametes, pop.mutations);
        }

        inline void
        update(const fwdpy11::SlocusPop& pop)
        {
            gv2w->update(pop);
            noise_fxn->update(pop);
        }
    };
} // namespace fwdpy11
#endif
