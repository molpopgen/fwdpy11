#include "diploid_pop_fitness.hpp"

void
calculate_diploid_fitness(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    const std::vector<fwdpy11::DiploidPopulationGeneticValue *>
        &gvalue_pointers,
    const std::vector<std::size_t> &deme_to_gvalue_map,
    std::vector<fwdpy11::DiploidMetadata> &offspring_metadata,
    std::vector<double> &new_diploid_gvalues,
    const bool update_genotype_matrix)
{
    // Calculate parental fitnesses
    double sum_parental_fitnesses = 0.0;
    if (update_genotype_matrix == true)
        {
            new_diploid_gvalues.resize(offspring_metadata.size()
                                       * gvalue_pointers[0]->total_dim);
        }
    auto gvoffset = new_diploid_gvalues.data();
    for (std::size_t i = 0; i < pop.diploids.size();
         ++i, gvoffset += gvalue_pointers[0]->total_dim)
        {
            auto idx = deme_to_gvalue_map[offspring_metadata[i].deme];
            gvalue_pointers[idx]->operator()(rng, i, pop,
                                             offspring_metadata[i]);
            if (update_genotype_matrix == true)
                {
                    std::copy(begin(gvalue_pointers[idx]->gvalues),
                              end(gvalue_pointers[idx]->gvalues), gvoffset);
                }
            sum_parental_fitnesses += offspring_metadata[i].w;
        }
    // If the sum of parental fitnesses is not finite,
    // then the genetic value calculator returned a non-finite value/
    // Unfortunately, gsl_ran_discrete_preproc allows such values through
    // without raising an error, so we have to check things here.
    if (!std::isfinite(sum_parental_fitnesses))
        {
            throw std::runtime_error("non-finite fitnesses encountered");
        }
}

fwdpp::gsl_ran_discrete_t_ptr
calculate_diploid_fitness_genomes(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
    std::vector<fwdpy11::DiploidMetadata> &offspring_metadata)
{
    // Calculate parental fitnesses
    std::vector<double> parental_fitnesses(pop.diploids.size());
    double sum_parental_fitnesses = 0.0;
    for (std::size_t i = 0; i < pop.diploids.size(); ++i)
        {
            genetic_value_fxn(rng, i, pop, offspring_metadata[i]);
            parental_fitnesses[i] = offspring_metadata[i].w;
            sum_parental_fitnesses += parental_fitnesses[i];
        }
    // If the sum of parental fitnesses is not finite,
    // then the genetic value calculator returned a non-finite value/
    // Unfortunately, gsl_ran_discrete_preproc allows such values through
    // without raising an error, so we have to check things here.
    if (!std::isfinite(sum_parental_fitnesses))
        {
            throw std::runtime_error("non-finite fitnesses encountered");
        }

    auto rv = fwdpp::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(
        parental_fitnesses.size(), parental_fitnesses.data()));
    if (rv == nullptr)
        {
            // This is due to negative fitnesses
            throw std::runtime_error(
                "fitness lookup table could not be generated");
        }
    return rv;
}
