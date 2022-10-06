#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace
{
    class GBR : public fwdpy11::DiploidGeneticValue
    {
      private:
        mutable std::size_t deme;
        template <typename mutation_key_cont_t, typename mcont_t>
        inline double
        sum_haplotype_effect_sizes(const mutation_key_cont_t& keys,
                                   const mcont_t& mutations) const
        {
            double rv = 0.0;
            for (auto& k : keys)
                {
                    rv += mutations[k].s;
                }
            return rv;
        }

        template <typename mutation_key_cont_t, typename mcont_t>
        inline double
        sum_haplotype_effect_sizes(std::size_t deme, const mutation_key_cont_t& keys,
                                   const mcont_t& mutations) const
        {
            double rv = 0.0;
            for (auto& k : keys)
                {
                    rv += mutations[k].esizes[deme];
                }
            return rv;
        }

        std::function<double(const std::size_t diploid_index,
                             const fwdpy11::DiploidPopulation& pop)>
        generate_backend_function(std::size_t ndemes, std::size_t& deme)
        {
            if (ndemes == 1)
                {
                    return [this](const std::size_t diploid_index,
                                  const fwdpy11::DiploidPopulation& pop) {
                        double h1 = sum_haplotype_effect_sizes(
                            pop.haploid_genomes[pop.diploids[diploid_index].first]
                                .smutations,
                            pop.mutations);
                        double h2 = sum_haplotype_effect_sizes(
                            pop.haploid_genomes[pop.diploids[diploid_index].second]
                                .smutations,
                            pop.mutations);
                        return sqrt(h1 * h2);
                    };
                }
            return [&deme, this](const std::size_t diploid_index,
                                 const fwdpy11::DiploidPopulation& pop) {
                double h1 = sum_haplotype_effect_sizes(
                    deme,
                    pop.haploid_genomes[pop.diploids[diploid_index].first].smutations,
                    pop.mutations);
                double h2 = sum_haplotype_effect_sizes(
                    deme,
                    pop.haploid_genomes[pop.diploids[diploid_index].second].smutations,
                    pop.mutations);
                return sqrt(h1 * h2);
            };
        }

        const std::function<double(const std::size_t diploid_index,
                                   const fwdpy11::DiploidPopulation& pop)>
            f;

      public:
        GBR(std::size_t ndim, const fwdpy11::GeneticValueToFitnessMap* gv2w_,
            const fwdpy11::GeneticValueNoise* noise_)
            : fwdpy11::DiploidGeneticValue{ndim, gv2w_, noise_},
              deme(0), f{generate_backend_function(total_dim, deme)}
        {
        }

        double
        calculate_gvalue(const fwdpy11::DiploidGeneticValueData data) override
        {
            deme = data.pop.get()
                       .diploid_metadata[data.offspring_metadata.get().label]
                       .deme;
            gvalues[0] = f(data.offspring_metadata.get().label, data.pop.get());
            return gvalues[0];
        }

        void
        update(const fwdpy11::DiploidPopulation& /*pop*/) override
        {
        }
    };
}

void
init_GBR(py::module& m)
{
    py::class_<GBR, fwdpy11::DiploidGeneticValue>(m, "_ll_GBR")
        .def(py::init([](const fwdpy11::GeneticValueIsTrait* gvalue_to_fitness,
                         const fwdpy11::GeneticValueNoise* noise) {
                 return GBR(1, gvalue_to_fitness, noise);
             }),
             py::arg("gvalue_to_fitness"), py::arg("noise"));
}
