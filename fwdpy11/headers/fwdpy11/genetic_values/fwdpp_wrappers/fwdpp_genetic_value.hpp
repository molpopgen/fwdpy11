#ifndef FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP__GVALUE_HPP__
#define FWDPY11_GENETIC_VALUES_WRAPPERS_FWDPP__GVALUE_HPP__

#include <type_traits>
#include <functional>
#include "../DiploidGeneticValue.hpp"
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>

namespace fwdpy11
{
    template <typename single_deme_het_fxn, typename single_deme_hom_fxn,
              typename multi_deme_het_fxn, typename multi_deme_hom_fxn,
              typename pickle_fxn, int starting_value>
    class stateless_site_dependent_genetic_value_wrapper : public DiploidGeneticValue
    {
      private:
        struct single_deme_callback
        {
            double
            operator()(const stateless_site_dependent_genetic_value_wrapper* outer_this,
                       const std::size_t diploid_index, const DiploidPopulation& pop,
                       const DiploidMetadata& /*metadata*/) const
            {
                outer_this->gvalues[0] = outer_this->make_return_value(
                    outer_this->gv(pop.diploids[diploid_index], pop.haploid_genomes,
                                   pop.mutations, outer_this->single_deme_aa,
                                   outer_this->single_deme_Aa, starting_value));
                return outer_this->gvalues[0];
            }
        };

        struct multi_deme_callback
        {
            double
            operator()(const stateless_site_dependent_genetic_value_wrapper* outer_this,
                       const std::size_t diploid_index, const DiploidPopulation& pop,
                       const DiploidMetadata& metadata) const
            {
                std::size_t deme = metadata.deme;
                outer_this->gvalues[0] = outer_this->make_return_value(outer_this->gv(
                    pop.diploids[diploid_index], pop.haploid_genomes, pop.mutations,
                    [deme, outer_this](double& d, const Mutation& mut) {
                        return outer_this->multi_deme_aa(deme, d, mut);
                    },
                    [deme, outer_this](double& d, const Mutation& mut) {
                        return outer_this->multi_deme_Aa(deme, d, mut);
                    },
                    starting_value));
                return outer_this->gvalues[0];
            }
        };

        using callback_type = std::function<double(
            const stateless_site_dependent_genetic_value_wrapper*, const std::size_t,
            const DiploidPopulation&, const DiploidMetadata&)>;

        callback_type
        init_callback(std::size_t n)
        {
            if (n == 1)
                {
                    return single_deme_callback();
                }
            return multi_deme_callback();
        }

        fwdpp::site_dependent_genetic_value gv;
        double aa_scaling;
        decltype(single_deme_het_fxn()) single_deme_Aa;
        decltype(single_deme_hom_fxn(aa_scaling)) single_deme_aa;
        decltype(multi_deme_het_fxn()) multi_deme_Aa;
        decltype(multi_deme_hom_fxn(aa_scaling)) multi_deme_aa;
        std::function<double(double)> make_return_value;
        callback_type callback;
        bool isfitness;

      public:
        template <typename make_return_value_fxn>
        stateless_site_dependent_genetic_value_wrapper(std::size_t ndim, double scaling,
                                                       make_return_value_fxn&& mrv)
            : DiploidGeneticValue{ndim}, gv{}, aa_scaling(scaling),
              single_deme_Aa(single_deme_het_fxn()),
              single_deme_aa(single_deme_hom_fxn(aa_scaling)),
              multi_deme_Aa(multi_deme_het_fxn()),
              multi_deme_aa(multi_deme_hom_fxn(aa_scaling)),
              make_return_value(std::forward<make_return_value_fxn>(mrv)),
              callback(init_callback(ndim)), isfitness(true)
        {
        }

        // NOTE: the following two constructors ASSUME
        // that isfitness == false!!!

        template <typename make_return_value_fxn>
        stateless_site_dependent_genetic_value_wrapper(
            std::size_t ndim, double scaling, make_return_value_fxn&& mrv,
            const GeneticValueToFitnessMap& gv2w_)
            : DiploidGeneticValue{ndim, gv2w_}, gv{}, aa_scaling(scaling),
              single_deme_Aa(single_deme_het_fxn()),
              single_deme_aa(single_deme_hom_fxn(aa_scaling)),
              multi_deme_Aa(multi_deme_het_fxn()),
              multi_deme_aa(multi_deme_hom_fxn(aa_scaling)),
              make_return_value(std::forward<make_return_value_fxn>(mrv)),
              callback(init_callback(ndim)), isfitness(false)
        {
        }

        template <typename make_return_value_fxn>
        stateless_site_dependent_genetic_value_wrapper(
            std::size_t ndim, double scaling, make_return_value_fxn&& mrv,
            const GeneticValueToFitnessMap& gv2w_, const GeneticValueNoise& noise_)
            : DiploidGeneticValue{ndim, gv2w_, noise_}, gv{}, aa_scaling(scaling),
              single_deme_Aa(single_deme_het_fxn()),
              single_deme_aa(single_deme_hom_fxn(aa_scaling)),
              multi_deme_Aa(multi_deme_het_fxn()),
              multi_deme_aa(multi_deme_hom_fxn(aa_scaling)),
              make_return_value(std::forward<make_return_value_fxn>(mrv)),
              callback(init_callback(ndim)), isfitness(false)
        {
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const DiploidMetadata& metadata,
                         const DiploidPopulation& pop) const override
        {
            return callback(this, diploid_index, pop, metadata);
        }

        pybind11::object
        pickle() const override
        {
            return pickle_fxn()(*this);
        }

        double
        scaling() const
        {
            return aa_scaling;
        }

        bool
        is_fitness() const
        {
            return isfitness;
        }
    };
} // namespace fwdpy11
#endif

