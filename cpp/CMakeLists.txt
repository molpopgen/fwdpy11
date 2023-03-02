set(FWDPP_TYPES_SOURCES fwdpp_types/init.cc ;
    fwdpp_types/MutationBase.cc
    fwdpp_types/HaploidGenome.cc
    fwdpp_types/DataMatrix.cc
    fwdpp_types/NULL_NODE.cc
    fwdpp_types/Node.cc
    fwdpp_types/Edge.cc
    fwdpp_types/MutationRecord.cc
    fwdpp_types/Site.cc
    fwdpp_types/NodeTable.cc
    fwdpp_types/EdgeTable.cc
    fwdpp_types/MutationTable.cc
    fwdpp_types/SiteTable.cc
    fwdpp_types/TableCollection.cc)


set(FWDPY11_TYPES_SOURCES fwdpy11_types/init.cc
    fwdpy11_types/Mutation.cc
    fwdpy11_types/MutationVector.cc
    fwdpy11_types/DiploidGenotype.cc
    fwdpy11_types/DiploidMetadata.cc
    fwdpy11_types/DiploidVector.cc
    fwdpy11_types/HaploidGenomeVector.cc
    fwdpy11_types/rng.cc
    fwdpy11_types/PopulationBase.cc
    fwdpy11_types/DiploidPopulation.cc
    fwdpy11_types/ts_from_tskit.cc
    fwdpy11_types/tsrecorders.cc
    fwdpy11_types/RecordNothing.cc
    fwdpy11_types/GeneticMapUnit.cc)

set(REGION_SOURCES regions/init.cc regions/Region.cc regions/Sregion.cc regions/GammaS.cc regions/ConstantS.cc
    regions/ExpS.cc regions/UniformS.cc regions/GaussianS.cc regions/MutationRegions.cc
    regions/RecombinationRegions.cc regions/MultivariateGaussianEffects.cc
    regions/mvDES.cc
    regions/LogNormalS.cc
    regions/DiscreteDESD.cc)

set(GSL_SOURCES gsl/init.cc
    gsl/gsl_random.cc)

set(GENETIC_VALUE_SOURCES genetic_values/init.cc
    genetic_values/DiploidGeneticValue.cc
    genetic_values/PyDiploidGeneticValue.cc
    genetic_values/Additive.cc
    genetic_values/Multiplicative.cc
    genetic_values/GBR.cc
    genetic_values/DiploidMultivariateEffectsStrictAdditive.cc
    genetic_values/dgvalue_pointer_vector.cc)

set(GENETIC_VALUE_TO_FITNESS_SOURCES
    genetic_value_to_fitness/init.cc
    genetic_value_to_fitness/GeneticValueToFitnessMap.cc
    genetic_value_to_fitness/GeneticValueIsFitness.cc
    genetic_value_to_fitness/GeneticValueIsTrait.cc
    genetic_value_to_fitness/GSSmo.cc
    genetic_value_to_fitness/MultivariateGSSmo.cc
    genetic_value_to_fitness/Optimum.cc
    genetic_value_to_fitness/PleiotropicOptima.cc)

set(GENETIC_VALUE_NOISE_SOURCES
    genetic_value_noise/init.cc
    genetic_value_noise/GeneticValueNoise.cc
    genetic_value_noise/NoNoise.cc
    genetic_value_noise/GaussianNoise.cc)

set(TS_SOURCES ts/init.cc ts/TreeIterator.cc ts/VariantIterator.cc
    ts/count_mutations.cc
    ts/simplify.cc
    ts/data_matrix_from_tables.cc
    ts/infinite_sites.cc
    ts/DataMatrixIterator.cc
    ts/node_traversal.cc)

set(EVOLVE_POPULATION_SOURCES evolve_population/init.cc
    evolve_population/_evolvets.cc
    evolve_population/no_stopping.cc)

set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/init.cc
    discrete_demography/exceptions.cc)

set (ARRAY_PROXY_SOURCES array_proxies/init.cc)

set (MUTATION_DOMINANCE_SOURCES mutation_dominance/init.cc
     mutation_dominance/MutationDominance.cc)

set (FUNCTION_SOURCES functions/init.cc
     functions/add_mutation.cc
     functions/_add_mutation.cc)

set (DEMES_SOURCES
    demes/init.cc
    demes/forward_graph.cc)

set(ALL_SOURCES ${FWDPP_TYPES_SOURCES}
    ${FWDPY11_TYPES_SOURCES}
    ${REGION_SOURCES}
    ${GENETIC_VALUE_NOISE_SOURCES}
    ${GENETIC_VALUE_TO_FITNESS_SOURCES}
    ${GENETIC_VALUE_SOURCES}
    ${TS_SOURCES}
    ${GSL_SOURCES}
    ${EVOLVE_POPULATION_SOURCES}
    ${DISCRETE_DEMOGRAPHY_SOURCES}
    ${ARRAY_PROXY_SOURCES}
    ${MUTATION_DOMINANCE_SOURCES}
    ${FUNCTION_SOURCES}
    ${DEMES_SOURCES})

# NOTE: unclear if I can assign this to a
# named target?
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(LTO_OPTIONS)
if(ENABLE_PROFILING OR DISABLE_LTO)
    set(LTO_OPTIONS NO_EXTRAS)
endif()
# These are the main modules
pybind11_add_module(_fwdpy11 MODULE ${LTO_OPTIONS} _fwdpy11.cc ${ALL_SOURCES})
add_dependencies(_fwdpy11 fwdpy11core)
set_target_properties(_fwdpy11 fwdpy11core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_directories(_fwdpy11 PRIVATE ${CMAKE_SOURCE_DIR}/fwdpy11)
set_target_properties(_fwdpy11 PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
# From https://github.com/pybind/cmake_example/issues/11
if (UNIX AND NOT APPLE)
    set_target_properties(_fwdpy11 PROPERTIES INSTALL_RPATH "$\{ORIGIN\}")
elseif(APPLE)
    set_target_properties(_fwdpy11 PROPERTIES INSTALL_RPATH "@loader_path")
endif()

if (ENABLE_PROFILING)
    set_target_properties(_fwdpy11 PROPERTIES CXX_VISIBILITY_PRESET "default")
endif()
target_link_libraries(_fwdpy11 PRIVATE GSL::gsl GSL::gslcblas fwdpy11core)
