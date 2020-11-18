#ifndef FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP
#define FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP

#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>

// Changed in 0.6.0 to return void, as sims w/tree
// sequences generate fitness lookups via DiscreteDemography
void
calculate_diploid_fitness(const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
                          std::vector<fwdpy11::DiploidGeneticValue *> &gvalue_pointers,
                          const std::vector<std::size_t> &deme_to_gvalue_map,
                          std::vector<fwdpy11::DiploidMetadata> &offspring_metadata,
                          std::vector<double> &new_diploid_gvalues,
                          const bool update_genotype_matrix);

#endif
