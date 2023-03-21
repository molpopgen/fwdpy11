set(MUTATION_DOMINANCE_SOURCES mutation_dominance/dependency_injection.cc)

set(DEMES_SOURCES
    demes/forward_graph.cc)

set(EVOLVE_DISCRETE_DEMES_SOURCES
    evolve_discrete_demes/cleanup_metadata.cc
    evolve_discrete_demes/diploid_pop_fitness.cc
    evolve_discrete_demes/evolvets.cc
    evolve_discrete_demes/index_and_count_mutations.cc
    evolve_discrete_demes/remove_extinct_genomes.cc
    evolve_discrete_demes/remove_extinct_mutations.cc
    evolve_discrete_demes/track_ancestral_counts.cc
    evolve_discrete_demes/track_mutation_counts.cc
    evolve_discrete_demes/runtime_checks.cc
    evolve_discrete_demes/util.cc
    evolve_discrete_demes/discrete_demography/simulation/pick_parents.cc
    evolve_discrete_demes/discrete_demography/simulation/validate_parental_state.cc)

set(GENETIC_MAP_SOURCES
    genetic_maps/regions.cc)

set(DIPLOID_POPULATION_SOURCES
    diploid_population/set_mutations.cc)

set(ALL_SOURCES
    ${MUTATION_DOMINANCE_SOURCES}
    ${DEMES_SOURCES}
    ${GENETIC_MAP_SOURCES}
    ${DIPLOID_POPULATION_SOURCES}
    ${EVOLVE_DISCRETE_DEMES_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
add_library(fwdpy11core SHARED ${ALL_SOURCES})
add_dependencies(fwdpy11core cargo-build_fp11rust header)
target_link_libraries(fwdpy11core LINK_PRIVATE ${CMAKE_BINARY_DIR}/rust/libfp11rust.a)
target_link_libraries(fwdpy11core PRIVATE GSL::gsl GSL::gslcblas)
