find_package(GSL REQUIRED)
option(USE_WEFFCPP "Use -Weffc++ during compilation" ON)
include_directories(BEFORE ${fwdpy11_SOURCE_DIR}/fwdpy11/headers ${fwdpy11_SOURCE_DIR}/fwdpy11/headers/fwdpp)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})

if (USE_WEFFCPP)
    add_compile_options(-Weffc++)
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

set(REGION_SOURCES src/regions/init.cc src/regions/Region.cc src/regions/Sregion.cc src/regions/GammaS.cc src/regions/ConstantS.cc
    src/regions/ExpS.cc src/regions/UniformS.cc src/regions/GaussianS.cc src/regions/MutationRegions.cc
    src/regions/RecombinationRegions.cc src/regions/MultivariateGaussianEffects.cc
    src/regions/GeneticMapUnit.cc src/regions/PoissonInterval.cc 
    src/regions/BinomialPoint.cc
    src/regions/PoissonPoint.cc
    src/regions/FixedCrossovers.cc)

set(GENETIC_VALUE_SOURCES src/genetic_values/init.cc src/genetic_values/SlocusMultivariateEffectsStrictAdditive.cc
    src/genetic_values/MultivariateGeneticValueToFitnessMap.cc
    src/genetic_values/SlocusMultivariateGeneticValueWithMapping.cc src/genetic_values/MultivariateGSS.cc
    src/genetic_values/MultivariateGSSmo.cc)

set(TS_SOURCES src/ts/init.cc src/ts/tree_iterator.cc src/ts/VariantIterator.cc
    src/ts/count_mutations.cc
    src/ts/simplify.cc
    src/ts/data_matrix_from_tables.cc)

set(TSEVOLUTION_SOURCES src/tsevolution/tsevolution.cc
    src/tsevolution/slocuspop.cc
    src/tsevolution/util.cc
    src/tsevolution/cleanup_metadata.cc
    src/tsevolution/genetic_value_common.cc
    src/tsevolution/slocus_fitness.cc
    src/tsevolution/index_and_count_mutations.cc
    src/tsevolution/track_mutation_counts.cc
    src/tsevolution/no_stopping.cc
    src/tsevolution/remove_extinct_mutations.cc)

# These are the main modules
pybind11_add_module(_fwdpy11 MODULE src/_fwdpy11.cc ${REGION_SOURCES})
target_link_libraries(_fwdpy11 PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_init MODULE src/_init.cc)
target_link_libraries(_init PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(fwdpp_types MODULE src/fwdpp_types.cc)
target_link_libraries(fwdpp_types PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_opaque_gametes MODULE src/_opaque_gametes.cc)
target_link_libraries(_opaque_gametes PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_opaque_mutations MODULE src/_opaque_mutations.cc)
target_link_libraries(_opaque_mutations PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_opaque_diploids MODULE src/_opaque_diploids.cc)
target_link_libraries(_opaque_diploids PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_Population MODULE src/_Population.cc)
# target_link_libraries(_Population PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_Populations MODULE src/_Populations.cc)
target_link_libraries(_Populations PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(fwdpy11_types MODULE src/fwdpy11_types.cc)
target_link_libraries(fwdpy11_types PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(_gslrng MODULE src/_gslrng.cc)
target_link_libraries(_gslrng PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(sampling MODULE src/fwdpy11_sampling.cc)
target_link_libraries(sampling PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(genetic_values MODULE src/genetic_values.cc ${GENETIC_VALUE_SOURCES})
target_link_libraries(genetic_values PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(genetic_value_noise MODULE src/genetic_value_noise.cc)
target_link_libraries(genetic_value_noise PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(wright_fisher_slocus MODULE src/wright_fisher_slocus.cc)
target_link_libraries(wright_fisher_slocus PRIVATE GSL::gsl GSL::gslcblas)
pybind11_add_module(gsl_random MODULE src/gsl_random.cc)
target_link_libraries(gsl_random PRIVATE GSL::gsl GSL::gslcblas)
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
