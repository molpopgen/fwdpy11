#include <functional>
#include <cmath>
#include <stdexcept>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>
#include <fwdpy11/evolvets/sample_recorder_types.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
#include <fwdpy11/samplers.hpp>

void evolve_with_tree_sequences(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    fwdpy11::SampleRecorder &sr, const unsigned simplification_interval,
    fwdpy11::discrete_demography::DiscreteDemography &demography,
    const std::uint32_t simlen, const double mu_neutral, const double mu_selected,
    const fwdpy11::MutationRegions &mmodel, const fwdpy11::GeneticMap &rmodel,
    // NOTE: gvalue_pointers is a change in 0.6.0,
    // and the object holds non-const bare pointers
    // to objects owned by Python.
    fwdpy11::dgvalue_pointer_vector_ &gvalue_pointers,
    fwdpy11::DiploidPopulation_sample_recorder recorder,
    std::function<bool(const fwdpy11::DiploidPopulation &, const bool)>
        &stopping_criteron,
    // NOTE: this is the complement of what a user will input, which is "prune_selected"
    const bool preserve_selected_fixations, const bool suppress_edge_table_indexing,
    bool record_genotype_matrix, const bool track_mutation_counts_during_sim,
    const bool remove_extinct_mutations_at_finish,
    const bool reset_treeseqs_to_alive_nodes_after_simplification,
    const bool preserve_first_generation,
    const fwdpy11::DiploidPopulation_temporal_sampler &post_simplification_recorder);

