#ifndef FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP
#define FWDPY11_TSEVOLVE_SLOCUS_FITNESS_HPP

#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>

std::function<fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
    const fwdpy11::GSLrng_t &g, fwdpy11::SlocusPop &,
    const fwdpy11::SlocusPopGeneticValue &,
    std::vector<fwdpy11::DiploidMetadata> &, std::vector<double> &)>
wrap_calculate_fitness_SlocusPop(bool update_genotype_matrix);

#endif
