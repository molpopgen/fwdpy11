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
              int starting_value>
    class stateless_site_dependent_genetic_value_wrapper : public DiploidGeneticValue
    {
      private:
        struct single_deme_callback
        {
            const single_deme_het_fxn single_deme_Aa;
            const single_deme_hom_fxn single_deme_aa;

            explicit single_deme_callback(double aa_scaling)
                : single_deme_Aa(), single_deme_aa(aa_scaling)
            {
            }

            inline double
            operator()(const fwdpp::site_dependent_genetic_value& gv,
                       const std::size_t diploid_index,
                       const DiploidMetadata& /*metadata*/,
                       const DiploidPopulation& pop) const
            {
                return gv(pop.diploids[diploid_index], pop.haploid_genomes,
                          pop.mutations, single_deme_aa, single_deme_Aa, starting_value);
            }
        };

        struct multi_deme_callback
        {
            const multi_deme_het_fxn multi_deme_Aa;
            const multi_deme_hom_fxn multi_deme_aa;

            explicit multi_deme_callback(double aa_scaling)
                : multi_deme_Aa(), multi_deme_aa(aa_scaling)
            {
            }

            inline double
            operator()(const fwdpp::site_dependent_genetic_value& gv,
                       const std::size_t diploid_index, const DiploidMetadata& metadata,
                       const DiploidPopulation& pop) const
            {
                std::size_t deme = metadata.deme;
                return gv(
                    pop.diploids[diploid_index], pop.haploid_genomes, pop.mutations,
                    [deme, this](double& d, const Mutation& mut) {
                        return multi_deme_aa(deme, d, mut);
                    },
                    [deme, this](double& d, const Mutation& mut) {
                        return multi_deme_Aa(deme, d, mut);
                    },
                    starting_value);
            }
        };

        using callback_type = std::function<double(
            const fwdpp::site_dependent_genetic_value&, const std::size_t,
            const DiploidMetadata&, const DiploidPopulation&)>;

        using make_return_value_t = std::function<double(double)>;

        callback_type
        init_callback(std::size_t n, double aa_scaling)
        {
            if (n == 1)
                {
                    return single_deme_callback(aa_scaling);
                }
            return multi_deme_callback(aa_scaling);
        }

        fwdpp::site_dependent_genetic_value gv;
        double aa_scaling;
        make_return_value_t make_return_value;
        callback_type callback;
        bool isfitness;

      public:
        stateless_site_dependent_genetic_value_wrapper(
            std::size_t ndim, double scaling, make_return_value_t mrv,
            const GeneticValueToFitnessMap* gv2w_, const GeneticValueNoise* noise_)
            : DiploidGeneticValue{ndim, gv2w_, noise_}, gv{}, aa_scaling(scaling),
              make_return_value(std::move(mrv)),
              callback(init_callback(ndim, aa_scaling)), isfitness(gv2w->isfitness)
        {
        }

        double
        calculate_gvalue(const DiploidGeneticValueData data) override
        {
            gvalues[0] = make_return_value(
                callback(gv, data.offspring_metadata.get().label,
                         data.offspring_metadata.get(), data.pop.get()));
            return gvalues[0];
        }

        void
        update(const fwdpy11::DiploidPopulation& /*pop*/) override
        {
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

