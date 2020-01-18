#ifndef FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP
#define FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP

#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>

// Changed in 0.6.0 to return void, as sims w/tree
// sequences generate fitness lookups via DiscreteDemography
void calculate_diploid_fitness(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    const std::vector<fwdpy11::DiploidPopulationGeneticValue *>
        &gvalue_pointers,
    const std::vector<std::size_t> &deme_to_gvalue_map,
    std::vector<fwdpy11::DiploidMetadata> &new_metadata,
    std::vector<double> &new_diploid_gvalues,
    const bool update_genotype_matrix);

// This overload was added in 0.6.0 as a temporary
// hack b/c sims with tree sequences generate fitness
// lookups via DiscreteDemography
fwdpp::gsl_ran_discrete_t_ptr calculate_diploid_fitness_genomes(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    const fwdpy11::DiploidPopulationGeneticValue &genetic_value_fxn,
    std::vector<fwdpy11::DiploidMetadata> &new_metadata);

#endif
