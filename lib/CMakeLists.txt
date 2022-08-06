set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/SetDemeSize.cc)

set(MUTATION_DOMINANCE_SOURCES mutation_dominance/dependency_injection.cc)

set(ALL_SOURCES ${DISCRETE_DEMOGRAPHY_SOURCES}
    ${MUTATION_DOMINANCE_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
add_library(fwdpy11core SHARED ${ALL_SOURCES})
target_link_libraries(fwdpy11core PRIVATE GSL::gsl GSL::gslcblas)
