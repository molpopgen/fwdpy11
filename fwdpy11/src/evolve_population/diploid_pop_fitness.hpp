#ifndef FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP
#define FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP

#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>

std::function<fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
    const fwdpy11::GSLrng_t &g, fwdpy11::DiploidPopulation &,
    const fwdpy11::DiploidPopulationGeneticValue &,
    std::vector<fwdpy11::DiploidMetadata> &, std::vector<double> &)>
wrap_calculate_fitness_DiploidPopulation(bool update_genotype_matrix);

#endif
