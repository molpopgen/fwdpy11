set(FWDPP_TYPES_SOURCES src/fwdpp_types/init.cc 
    src/fwdpp_types/MutationBase.cc
    src/fwdpp_types/HaploidGenome.cc
    src/fwdpp_types/DataMatrix.cc
    src/fwdpp_types/NULL_NODE.cc
    src/fwdpp_types/Node.cc
    src/fwdpp_types/Edge.cc
    src/fwdpp_types/MutationRecord.cc
    src/fwdpp_types/Site.cc
    src/fwdpp_types/IndexedEdge.cc
    src/fwdpp_types/NodeTable.cc
    src/fwdpp_types/EdgeTable.cc
    src/fwdpp_types/MutationTable.cc
    src/fwdpp_types/SiteTable.cc
    src/fwdpp_types/TableCollection.cc
    src/fwdpp_types/GeneticMapUnit.cc 
    src/fwdpp_types/PoissonInterval.cc 
    src/fwdpp_types/BinomialPoint.cc
    src/fwdpp_types/PoissonPoint.cc
    src/fwdpp_types/FixedCrossovers.cc
    src/fwdpp_types/BinomialInterval.cc)


set(FWDPP_FUNCTIONS_SOURCES src/fwdpp_functions/init.cc
    src/fwdpp_functions/data_matrix_creation.cc)

set(FWDPY11_TYPES_SOURCES src/fwdpy11_types/init.cc
    src/fwdpy11_types/Mutation.cc
    src/fwdpy11_types/MutationVector.cc
    src/fwdpy11_types/DiploidGenotype.cc
    src/fwdpy11_types/DiploidMetadata.cc
    src/fwdpy11_types/DiploidVector.cc
    src/fwdpy11_types/HaploidGenomeVector.cc
    src/fwdpy11_types/rng.cc
    src/fwdpy11_types/PopulationBase.cc
    src/fwdpy11_types/DiploidPopulation.cc
    src/fwdpy11_types/ts_from_tskit.cc
    src/fwdpy11_types/tsrecorders.cc
    src/fwdpy11_types/RecordNothing.cc)

set(FWDPY11_FUNCTIONS_SOURCES src/fwdpy11_functions/init.cc
    src/fwdpy11_functions/change_effect_size.cc
    src/fwdpy11_functions/sort_gamete_keys.cc)

set(REGION_SOURCES src/regions/init.cc src/regions/Region.cc src/regions/Sregion.cc src/regions/GammaS.cc src/regions/ConstantS.cc
    src/regions/ExpS.cc src/regions/UniformS.cc src/regions/GaussianS.cc src/regions/MutationRegions.cc
    src/regions/RecombinationRegions.cc src/regions/MultivariateGaussianEffects.cc)

set(GSL_SOURCES src/gsl/init.cc
    src/gsl/gsl_random.cc)

set(GENETIC_VALUE_SOURCES src/genetic_values/init.cc
    src/genetic_values/GeneticValue.cc
    src/genetic_values/GeneticValueWithMapping.cc
    src/genetic_values/Additive.cc
    src/genetic_values/Multiplicative.cc
    src/genetic_values/GBR.cc
    src/genetic_values/DiploidMultivariateEffectsStrictAdditive.cc
    src/genetic_values/DiploidMultivariateGeneticValueWithMapping.cc
    src/genetic_values/dgvalue_pointer_vector.cc)

set(GENETIC_VALUE_TO_FITNESS_SOURCES
    src/genetic_value_to_fitness/init.cc
    src/genetic_value_to_fitness/GeneticValueToFitnessMap.cc
    src/genetic_value_to_fitness/GeneticValueIsFitness.cc
    src/genetic_value_to_fitness/GeneticValueIsTrait.cc
    src/genetic_value_to_fitness/GSS.cc
    src/genetic_value_to_fitness/GSSmo.cc
    src/genetic_value_to_fitness/MultivariateGeneticValueToFitnessMap.cc
    src/genetic_value_to_fitness/MultivariateGSS.cc
    src/genetic_value_to_fitness/MultivariateGSSmo.cc)

set(GENETIC_VALUE_NOISE_SOURCES
    src/genetic_value_noise/init.cc
    src/genetic_value_noise/GeneticValueNoise.cc
    src/genetic_value_noise/NoNoise.cc
    src/genetic_value_noise/GaussianNoise.cc)

set(TS_SOURCES src/ts/init.cc src/ts/TreeIterator.cc src/ts/VariantIterator.cc
    src/ts/count_mutations.cc
    src/ts/simplify.cc
    src/ts/data_matrix_from_tables.cc
    src/ts/infinite_sites.cc
    src/ts/DataMatrixIterator.cc
    src/ts/node_traversal.cc)

set(EVOLVE_POPULATION_SOURCES src/evolve_population/init.cc
    src/evolve_population/with_tree_sequences.cc
    src/evolve_population/no_tree_sequences.cc
    src/evolve_population/util.cc
    src/evolve_population/cleanup_metadata.cc
    src/evolve_population/diploid_pop_fitness.cc
    src/evolve_population/index_and_count_mutations.cc
    src/evolve_population/track_mutation_counts.cc
    src/evolve_population/no_stopping.cc
    src/evolve_population/remove_extinct_mutations.cc
    src/evolve_population/track_ancestral_counts.cc
    src/evolve_population/remove_extinct_genomes.cc)

set(DISCRETE_DEMOGRAPHY_SOURCES src/discrete_demography/init.cc
    src/discrete_demography/MigrationMatrix.cc
    src/discrete_demography/MassMigration.cc
    src/discrete_demography/SetDemeSize.cc
    src/discrete_demography/SetExponentialGrowth.cc
    src/discrete_demography/SetSelfingRate.cc
    src/discrete_demography/SetMigrationRates.cc
    src/discrete_demography/DiscreteDemography.cc
    src/discrete_demography/exceptions.cc)

set(ALL_SOURCES ${FWDPP_TYPES_SOURCES}
    ${FWDPP_FUNCTIONS_SOURCES}
    ${FWDPY11_TYPES_SOURCES}
    ${FWDPY11_FUNCTIONS_SOURCES}
    ${REGION_SOURCES}
    ${GENETIC_VALUE_NOISE_SOURCES}
    ${GENETIC_VALUE_TO_FITNESS_SOURCES}
    ${GENETIC_VALUE_SOURCES}
    ${TS_SOURCES}
    ${GSL_SOURCES}
    ${EVOLVE_POPULATION_SOURCES}
    ${DISCRETE_DEMOGRAPHY_SOURCES})

set(LTO_OPTIONS)
if(ENABLE_PROFILING OR DISABLE_LTO)
    set(LTO_OPTIONS NO_EXTRAS)
endif()
# These are the main modules
pybind11_add_module(_fwdpy11 MODULE ${LTO_OPTIONS} src/_fwdpy11.cc ${ALL_SOURCES})

if (ENABLE_PROFILING)
    set_target_properties(_fwdpy11 PROPERTIES CXX_VISIBILITY_PRESET "default")
endif()
target_link_libraries(_fwdpy11 PRIVATE GSL::gsl GSL::gslcblas)
