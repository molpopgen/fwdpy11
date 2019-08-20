cmake_minimum_required(VERSION 2.8.12)
project(gavlue_recorder)

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.2.3')
    message(FATAL_ERROR "pybind11 version must be >= '2.2.3'")
endif()

execute_process(COMMAND python3 -m fwdpy11 --fwdpy11_headers OUTPUT_VARIABLE FP11HEADERS)
execute_process(COMMAND python3 -m fwdpy11 --fwdpp_headers OUTPUT_VARIABLE FWDPPHEADERS)

find_package(GSL REQUIRED)
include_directories(BEFORE ${FP11HEADERS} ${FWDPPHEADERS})
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})


set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

pybind11_add_module(gvalue_recorder MODULE gvalue_recorder.cc)
target_link_libraries(gvalue_recorder PRIVATE GSL::gsl GSL::gslcblas)
