#ifndef FWDPY11_TSEVOLVE_MLOCUS_FITNESS_HPP
#define FWDPY11_TSEVOLVE_MLOCUS_FITNESS_HPP

#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/genetic_values/MlocusPopGeneticValue.hpp>

std::function<fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
    const fwdpy11::GSLrng_t &g, fwdpy11::MlocusPop &,
    const fwdpy11::MlocusPopGeneticValue &,
    std::vector<fwdpy11::DiploidMetadata> &, std::vector<double> &)>
wrap_calculate_fitness_MlocusPop(bool update_genotype_matrix);

#endif

