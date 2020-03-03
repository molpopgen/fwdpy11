#ifndef FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP__GVALUE_HPP__
#define FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP__GVALUE_HPP__

#include <type_traits>
#include <functional>
#include "../DiploidGeneticValue.hpp"
#include "../noise.hpp"

namespace fwdpy11
{
    template <typename fwdppT, typename pickleFunction>
    struct fwdpp_genetic_value : public fwdpy11::DiploidGeneticValue
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
        explicit fwdpp_genetic_value(forwarded_fwdppT&& gv_)
            : DiploidGeneticValue{ 1 }, gv{ std::forward<forwarded_fwdppT>(
                                            gv_) },
              pickle_fxn(pickleFunction{})
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_genetic_value(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_)
            : DiploidGeneticValue{ 1, gv2w_ },
              gv{ std::forward<forwarded_fwdppT>(gv_) },
              pickle_fxn(pickleFunction())
        {
        }

        template <typename forwarded_fwdppT>
        fwdpp_genetic_value(forwarded_fwdppT&& gv_,
                            const GeneticValueToFitnessMap& gv2w_,
                            const GeneticValueNoise& noise_)
            : DiploidGeneticValue{ 1, gv2w_, noise_ },
              gv{ std::forward<forwarded_fwdppT>(gv_)

              },
              pickle_fxn(pickleFunction())

        {
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const fwdpy11::DiploidPopulation& pop) const override
        {
            gvalues[0] = gv(pop.diploids[diploid_index], pop.haploid_genomes,
                            pop.mutations);
            return gvalues[0];
        }

        virtual pybind11::object
        pickle() const
        {
            return pickle_fxn(gv);
        }
    };
} // namespace fwdpy11
#endif
