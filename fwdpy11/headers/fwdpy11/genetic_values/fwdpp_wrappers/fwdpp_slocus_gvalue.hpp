#ifndef FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_SLOCUS_GVALUE_HPP__
#define FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP_SLOCUS_GVALUE_HPP__

#include <type_traits>
#include <functional>
#include "../SlocusPopGeneticValueWithMapping.hpp"
#include "../noise.hpp"

namespace fwdpy11
{
    template <typename fwdppT, typename pickleFunction>
    struct fwdpp_slocus_gvalue
        : public fwdpy11::SlocusPopGeneticValueWithMapping
    {
        using gvalue_map_ptr
            = std::unique_ptr<fwdpy11::GeneticValueToFitnessMap>;
        const fwdppT gv;
        const pickleFunction pickle_fxn;
        static_assert(
            std::is_convertible<pickleFunction, std::function<pybind11::object(
                                                    const fwdppT&)>>::value,
            "pickling function must be convertible to "
            "std::function<pybind11::object(const fwdppT*)>");

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_)
            : SlocusPopGeneticValueWithMapping{ GeneticValueIsFitness() },
              gv{ std::forward<forwarded_fwdppT>(gv_) }, pickle_fxn{
                  pickleFunction{}
              }
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_)
            : SlocusPopGeneticValueWithMapping{ gv2w_ },
              gv{ std::forward<forwarded_fwdppT>(gv_) }, pickle_fxn{
                  pickleFunction()
              }
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_slocus_gvalue(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : SlocusPopGeneticValueWithMapping{ gv2w_, noise_ },
              gv{ std::forward<forwarded_fwdppT>(gv_)

              },
              pickle_fxn{ pickleFunction() }

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

        virtual pybind11::object
        pickle() const
        {
            return pickle_fxn(gv);
        }
    };
} // namespace fwdpy11
#endif
