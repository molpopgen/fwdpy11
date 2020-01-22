cmake_minimum_required(VERSION 2.8.12)
project(fwdpy11)

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.2.3')
    message(FATAL_ERROR "pybind11 version must be >= '2.2.3'")
endif()

find_package(GSL REQUIRED)
option(USE_WEFFCPP "Use -Weffc++ during compilation" OFF)
option(ENABLE_PROFILING "Compile to enable code profiling" OFF)
option(BUILD_UNIT_TESTS "Build C++ modules for unit tests" ON)
option(DISABLE_LTO "Disable link-time optimization (LTO)" OFF)
include_directories(BEFORE ${fwdpy11_SOURCE_DIR}/fwdpy11/headers ${fwdpy11_SOURCE_DIR}/fwdpy11/headers/fwdpp)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})

if (USE_WEFFCPP)
    add_compile_options(-Weffc++)
endif()

if (ENABLE_PROFILING)
    add_definitions(-DPYBIND11_NAMESPACE=pybind11)
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
add_subdirectory(fwdpy11)

if(BUILD_UNIT_TESTS)
    add_subdirectory(tests)
endif(BUILD_UNIT_TESTS)
