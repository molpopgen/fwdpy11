SET(CPPTEST_SOURCES
    fwdpy11_cpp_tests.cc
    test_discrete_demography_stable_sorting_of_events.cc
    test_lowlevel_demographic_events.cc
    test_lowlevel_discrete_demography_bad_models.cc
    test_lowlevel_discrete_demography_objects.cc
    discrete_demography_roundtrips.cc
    discrete_demography_util.cc
    test_MutationDominance.cc
    test_gss_classes.cc
    test_effect_sizes_of_zero.cc
    test_Sregion_from_mvnorm.cc
    test_forward_demes_graph.cc
    forward_demes_graph_fixtures.cc
)

# This code is a WIP refactor
# to use the ForardGraph type that
# is backed by the demes-forward rust
# crate.
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
)

add_executable(fwdpy11_cpp_tests ${CPPTEST_SOURCES} ${EVOLVE_DISCRETE_DEMES_SOURCES})
add_dependencies(fwdpy11_cpp_tests fwdpy11core header)
target_link_directories(fwdpy11_cpp_tests PRIVATE ${CMAKE_SOURCE_DIR}/fwdpy11)
target_link_libraries(fwdpy11_cpp_tests PRIVATE boost_unit_test_framework fwdpy11core GSL::gsl GSL::gslcblas)
if (NOT APPLE)
    target_link_options(fwdpy11_cpp_tests BEFORE PUBLIC LINKER:--no-as-needed -ldl)
endif()
set_target_properties(fwdpy11_cpp_tests PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_TEST_DYN_LINK")

add_test(NAME fwdpy11_cpp_tests
    COMMAND fwdpy11_cpp_tests)
