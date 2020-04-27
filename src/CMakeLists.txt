set(FWDPP_TYPES_SOURCES fwdpp_types/init.cc 
    fwdpp_types/MutationBase.cc
    fwdpp_types/HaploidGenome.cc
    fwdpp_types/DataMatrix.cc
    fwdpp_types/NULL_NODE.cc
    fwdpp_types/Node.cc
    fwdpp_types/Edge.cc
    fwdpp_types/MutationRecord.cc
    fwdpp_types/Site.cc
    fwdpp_types/IndexedEdge.cc
    fwdpp_types/NodeTable.cc
    fwdpp_types/EdgeTable.cc
    fwdpp_types/MutationTable.cc
    fwdpp_types/SiteTable.cc
    fwdpp_types/TableCollection.cc
    fwdpp_types/GeneticMapUnit.cc 
    fwdpp_types/PoissonInterval.cc 
    fwdpp_types/BinomialPoint.cc
    fwdpp_types/PoissonPoint.cc
    fwdpp_types/FixedCrossovers.cc
    fwdpp_types/BinomialInterval.cc)


set(FWDPP_FUNCTIONS_SOURCES fwdpp_functions/init.cc
    fwdpp_functions/data_matrix_creation.cc)

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
    fwdpy11_types/RecordNothing.cc)

set(FWDPY11_FUNCTIONS_SOURCES fwdpy11_functions/init.cc
    fwdpy11_functions/change_effect_size.cc
    fwdpy11_functions/sort_gamete_keys.cc)

set(REGION_SOURCES regions/init.cc regions/Region.cc regions/Sregion.cc regions/GammaS.cc regions/ConstantS.cc
    regions/ExpS.cc regions/UniformS.cc regions/GaussianS.cc regions/MutationRegions.cc
    regions/RecombinationRegions.cc regions/MultivariateGaussianEffects.cc
    regions/mvDES.cc
    regions/LogNormalS.cc)

set(GSL_SOURCES gsl/init.cc
    gsl/gsl_random.cc)

set(GENETIC_VALUE_SOURCES genetic_values/init.cc
    genetic_values/DiploidGeneticValue.cc
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
    genetic_value_to_fitness/GSS.cc
    genetic_value_to_fitness/GSSmo.cc
    genetic_value_to_fitness/MultivariateGSS.cc
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
    evolve_population/with_tree_sequences.cc
    evolve_population/no_tree_sequences.cc
    evolve_population/util.cc
    evolve_population/cleanup_metadata.cc
    evolve_population/diploid_pop_fitness.cc
    evolve_population/index_and_count_mutations.cc
    evolve_population/track_mutation_counts.cc
    evolve_population/no_stopping.cc
    evolve_population/remove_extinct_mutations.cc
    evolve_population/track_ancestral_counts.cc
    evolve_population/remove_extinct_genomes.cc)

set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/init.cc
    discrete_demography/MigrationMatrix.cc
    discrete_demography/MassMigration.cc
    discrete_demography/SetDemeSize.cc
    discrete_demography/SetExponentialGrowth.cc
    discrete_demography/SetSelfingRate.cc
    discrete_demography/SetMigrationRates.cc
    discrete_demography/DiscreteDemography.cc
    discrete_demography/exceptions.cc)

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
pybind11_add_module(_fwdpy11 MODULE ${LTO_OPTIONS} _fwdpy11.cc ${ALL_SOURCES})

if (ENABLE_PROFILING)
    set_target_properties(_fwdpy11 PROPERTIES CXX_VISIBILITY_PRESET "default")
endif()
target_link_libraries(_fwdpy11 PRIVATE GSL::gsl GSL::gslcblas)
