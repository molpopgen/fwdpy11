set(DISCRETE_DEMOGRAPHY_SOURCES
    discrete_demography/SetDemeSize.cc)

set(MUTATION_DOMINANCE_SOURCES mutation_dominance/dependency_injection.cc)

set(DEMES_SOURCES
    demes/forward_graph.cc)

set(ALL_SOURCES ${DISCRETE_DEMOGRAPHY_SOURCES}
    ${MUTATION_DOMINANCE_SOURCES}
    ${DEMES_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
add_library(fwdpy11core SHARED ${ALL_SOURCES})
add_dependencies(fwdpy11core cargo-build_fp11rust header)
target_link_libraries(fwdpy11core LINK_PRIVATE ${CMAKE_BINARY_DIR}/rust/libfp11rust.a)
target_link_libraries(fwdpy11core PRIVATE GSL::gsl GSL::gslcblas)
