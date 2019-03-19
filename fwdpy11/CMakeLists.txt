find_package(GSL REQUIRED)
option(USE_WEFFCPP "Use -Weffc++ during compilation" ON)
include_directories(BEFORE ${fwdpy11_SOURCE_DIR}/fwdpy11/headers ${fwdpy11_SOURCE_DIR}/fwdpy11/headers/fwdpp)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})

if (USE_WEFFCPP)
    add_compile_options(-Weffc++)
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

set(FWDPP_TYPES_SOURCES src/fwdpp_types/init.cc 
    src/fwdpp_types/mutation_base.cc
    src/fwdpp_types/gamete.cc
    src/fwdpp_types/data_matrix.cc)

set(FWDPP_FUNCTIONS_SOURCES src/fwdpp_functions/init.cc
    src/fwdpp_functions/data_matrix_creation.cc)

set(FWDPY11_TYPES_SOURCES src/fwdpy11_types/init.cc
    src/fwdpy11_types/Mutation.cc
    src/fwdpy11_types/MutationVector.cc
    src/fwdpy11_types/DiploidGenotype.cc
    src/fwdpy11_types/DiploidMetadata.cc
    src/fwdpy11_types/DiploidVector.cc
    src/fwdpy11_types/GameteVector.cc
    src/fwdpy11_types/rng.cc)

set(REGION_SOURCES src/regions/init.cc src/regions/Region.cc src/regions/Sregion.cc src/regions/GammaS.cc src/regions/ConstantS.cc
    src/regions/ExpS.cc src/regions/UniformS.cc src/regions/GaussianS.cc src/regions/MutationRegions.cc
    src/regions/RecombinationRegions.cc src/regions/MultivariateGaussianEffects.cc
    src/regions/GeneticMapUnit.cc src/regions/PoissonInterval.cc 
    src/regions/BinomialPoint.cc
    src/regions/PoissonPoint.cc
    src/regions/FixedCrossovers.cc)

set(GSL_SOURCES src/gsl/init.cc
    src/gsl/gsl_random.cc)

set(GENETIC_VALUE_SOURCES src/genetic_values/init.cc src/genetic_values/DiploidMultivariateEffectsStrictAdditive.cc
    src/genetic_values/MultivariateGeneticValueToFitnessMap.cc
    src/genetic_values/DiploidMultivariateGeneticValueWithMapping.cc src/genetic_values/MultivariateGSS.cc
    src/genetic_values/MultivariateGSSmo.cc)

set(TS_SOURCES src/ts/init.cc src/ts/tree_iterator.cc src/ts/VariantIterator.cc
    src/ts/count_mutations.cc
    src/ts/simplify.cc
    src/ts/data_matrix_from_tables.cc)

set(TSEVOLUTION_SOURCES src/tsevolution/tsevolution.cc
    src/tsevolution/diploidpop.cc
    src/tsevolution/util.cc
    src/tsevolution/cleanup_metadata.cc
    src/tsevolution/genetic_value_common.cc
    src/tsevolution/diploid_pop_fitness.cc
    src/tsevolution/index_and_count_mutations.cc
    src/tsevolution/track_mutation_counts.cc
    src/tsevolution/no_stopping.cc
    src/tsevolution/remove_extinct_mutations.cc)

# These are the main modules
pybind11_add_module(_fwdpy11 MODULE src/_fwdpy11.cc ${FWDPP_TYPES_SOURCES}
    ${FWDPP_FUNCTIONS_SOURCES}
    ${FWDPY11_TYPES_SOURCES}
    ${REGION_SOURCES}
    ${GSL_SOURCES})
target_link_libraries(_fwdpy11 PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_Population MODULE src/_Population.cc)
# target_link_libraries(_Population PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_Populations MODULE src/_Populations.cc)
target_link_libraries(_Populations PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(genetic_values MODULE src/genetic_values.cc ${GENETIC_VALUE_SOURCES})
target_link_libraries(genetic_values PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(genetic_value_noise MODULE src/genetic_value_noise.cc)
target_link_libraries(genetic_value_noise PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_wright_fisher MODULE src/wright_fisher.cc)
target_link_libraries(_wright_fisher PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(util MODULE src/fwdpy11_util.cc)
target_link_libraries(util PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(ts MODULE src/ts.cc ${TS_SOURCES})
target_link_libraries(ts PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(tsrecorders MODULE src/tsrecorders.cc)
target_link_libraries(tsrecorders PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(ts_from_tskit MODULE src/ts_from_tskit.cc)
target_link_libraries(ts_from_tskit PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_tsevolution MODULE ${TSEVOLUTION_SOURCES})
target_link_libraries(_tsevolution PRIVATE GSL::gsl GSL::gslcblas)
