SET(CPP_NEUTRAL_BENCHMARK_SOURCES 
    cpp_neutral_benchmark.cc                                       
    ../fwdpy11/src/evolve_population/remove_extinct_genomes.cc
    ../fwdpy11/src/evolve_population/remove_extinct_mutations.cc
    ../fwdpy11/src/evolve_population/diploid_pop_fitness.cc
    ../fwdpy11/src/evolve_population/runtime_checks.cc
    ../fwdpy11/src/evolve_population/track_mutation_counts.cc
    ../fwdpy11/src/evolve_population/track_ancestral_counts.cc
    ../fwdpy11/src/evolve_population/index_and_count_mutations.cc
    ../fwdpy11/src/evolve_population/util.cc
    ../fwdpy11/src/evolve_population/evolvets.cc
    ../fwdpy11/src/evolve_population/cleanup_metadata.cc)

add_executable(cpp_neutral_benchmark ${CPP_NEUTRAL_BENCHMARK_SOURCES})
add_dependencies(cpp_neutral_benchmark fwdpy11core header)
target_link_directories(cpp_neutral_benchmark PRIVATE ${CMAKE_SOURCE_DIR}/fwdpy11)
target_link_libraries(cpp_neutral_benchmark PRIVATE fwdpy11core boost_program_options GSL::gsl GSL::gslcblas) 
if (NOT APPLE)
    target_link_options(cpp_neutral_benchmark BEFORE PUBLIC LINKER:--no-as-needed -ldl)
endif()

target_include_directories(cpp_neutral_benchmark BEFORE PUBLIC ${CMAKE_SOURCE_DIR}/fwdpy11/headers ${CMAKE_SOURCE_DIR}/fwdpy11/headers/fwdpp ${CMAKE_SOURCE_DIR}/fwdpy11/src/evolve_population)

