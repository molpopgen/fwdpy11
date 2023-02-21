#include "diploid_pop_fitness.hpp"

void
calculate_diploid_fitness(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    std::vector<fwdpy11::DiploidGeneticValue *> &gvalue_pointers,
    const std::vector<std::size_t> &deme_to_gvalue_map,
    std::vector<fwdpy11::DiploidMetadata> &offspring_metadata,
    std::vector<double> &new_diploid_gvalues, const bool update_genotype_matrix)
{
    // Calculate parental fitnesses
    double sum_parental_fitnesses = 0.0;
    new_diploid_gvalues.clear();
    for (std::size_t i = 0; i < offspring_metadata.size(); ++i)
        {
            auto idx = deme_to_gvalue_map[offspring_metadata[i].deme];
            //gvalue_pointers[idx]->operator()(rng, offspring_metadata[i].label, pop,
            //                                 offspring_metadata[i]);
            gvalue_pointers[idx]->operator()(fwdpy11::DiploidGeneticValueData(
                rng, pop, pop.diploid_metadata[offspring_metadata[i].parents[0]],
                pop.diploid_metadata[offspring_metadata[i].parents[1]], i,
                offspring_metadata[i]));
            if (update_genotype_matrix == true)
                {
                    new_diploid_gvalues.insert(end(new_diploid_gvalues),
                                               begin(gvalue_pointers[idx]->gvalues),
                                               end(gvalue_pointers[idx]->gvalues));
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
