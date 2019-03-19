#include "diploid_pop_fitness.hpp"
#include "genetic_value_common.hpp"

template <typename update_genotype_matrix>
fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitness_details(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
    std::vector<fwdpy11::DiploidMetadata> &new_metadata,
    std::vector<double> &new_diploid_gvalues, const update_genotype_matrix um)
{
    // Calculate parental fitnesses
    std::vector<double> parental_fitnesses(pop.diploids.size());
    double sum_parental_fitnesses = 0.0;
    new_metadata.resize(pop.N);
    resize_genotype_matrix(new_diploid_gvalues,
                           pop.N * genetic_value_fxn.total_dim, um);
    auto gvoffset = new_diploid_gvalues.data();
    for (std::size_t i = 0; i < pop.diploids.size();
         ++i, gvoffset += genetic_value_fxn.total_dim)
        {
            new_metadata[i] = pop.diploid_metadata[i];
            genetic_value_fxn(rng, i, pop, new_metadata[i]);
            copy_genetic_values(gvoffset, genetic_value_fxn.gvalues, um);
            parental_fitnesses[i] = new_metadata[i].w;
            sum_parental_fitnesses += parental_fitnesses[i];
        }
    pop.diploid_metadata.swap(new_metadata);
    pop.genetic_value_matrix.swap(new_diploid_gvalues);
    // If the sum of parental fitnesses is not finite,
    // then the genetic value calculator returned a non-finite value/
    // Unfortunately, gsl_ran_discrete_preproc allows such values through
    // without raising an error, so we have to check things here.
    if (!std::isfinite(sum_parental_fitnesses))
        {
            throw std::runtime_error("non-finite fitnesses encountered");
        }

    auto rv = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(parental_fitnesses.size(),
                                 parental_fitnesses.data()));
    if (rv == nullptr)
        {
            // This is due to negative fitnesses
            throw std::runtime_error(
                "fitness lookup table could not be generated");
        }
    return rv;
}

std::function<fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
    const fwdpy11::GSLrng_t &g, fwdpy11::DiploidPopulation &,
    const fwdpy11::DiploidPopulationGeneticValue &,
    std::vector<fwdpy11::DiploidMetadata> &, std::vector<double> &)>
wrap_calculate_fitness_DiploidPopulation(bool update_genotype_matrix)
{
    if (update_genotype_matrix)
        {
            return [](const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
                      const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
                      std::vector<fwdpy11::DiploidMetadata> &new_metadata,
                      std::vector<double> &new_diploid_gvalues) {
                return calculate_fitness_details(
                    rng, pop, genetic_value_fxn, new_metadata,
                    new_diploid_gvalues, std::true_type());
            };
        }
    return [](const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
              const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
              std::vector<fwdpy11::DiploidMetadata> &new_metadata,
              std::vector<double> &new_diploid_gvalues) {
        return calculate_fitness_details(rng, pop, genetic_value_fxn,
                                         new_metadata, new_diploid_gvalues,
                                         std::false_type());
    };
}

