set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/SetDemeSize.cc)

set(ALL_SOURCES ${DISCRETE_DEMOGRAPHY_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
add_library(fwdpy11core STATIC ${ALL_SOURCES})
