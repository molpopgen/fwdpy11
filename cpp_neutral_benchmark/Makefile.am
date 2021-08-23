noinst_PROGRAMS=cpp_neutral_benchmark

cpp_neutral_benchmark_SOURCES=cpp_neutral_benchmark.cc \
						  ../fwdpy11/src/evolve_population/evolvets.cc \
						  ../fwdpy11/src/evolve_population/util.cc \
						  ../fwdpy11/src/evolve_population/remove_extinct_genomes.cc \
						  ../fwdpy11/src/evolve_population/remove_extinct_mutations.cc \
						  ../fwdpy11/src/evolve_population/diploid_pop_fitness.cc \
						  ../fwdpy11/src/evolve_population/runtime_checks.cc \
						  ../fwdpy11/src/evolve_population/track_mutation_counts.cc \
						  ../fwdpy11/src/evolve_population/track_ancestral_counts.cc \
						  ../fwdpy11/src/evolve_population/index_and_count_mutations.cc \
						  ../fwdpy11/src/evolve_population/cleanup_metadata.cc

AM_CPPFLAGS=-I../fwdpy11/headers -I../fwdpy11/headers/fwdpp -I../fwdpy11/src/evolve_population
AM_CXXFLAGS=-W -Wall 

cpp_neutral_benchmark_LDADD=-lboost_program_options

AM_LIBS=

LIBS+=$(AM_LIBS)
