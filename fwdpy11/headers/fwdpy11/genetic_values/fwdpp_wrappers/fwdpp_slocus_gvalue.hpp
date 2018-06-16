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

        fwdpp_slocus_gvalue(const double);

        fwdpp_slocus_gvalue(const double scaling,
                                 const fwdpy11::GeneticValueIsTrait& g2w);

        fwdpp_slocus_gvalue(const double scaling,
                                 const fwdpy11::GeneticValueIsTrait& g2w,
                                 const fwdpy11::GeneticValueNoise& noise_fxn);

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
