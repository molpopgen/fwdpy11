cmake_minimum_required(VERSION 3.15)
project(fwdpy11)

set(CMAKE_CXX_EXTENSIONS OFF)
if(USECPP17)
    set(CMAKE_CXX_STANDARD 17) 
else()
    set(CMAKE_CXX_STANDARD 14) 
endif()
message(STATUS "Compiling with C++${CMAKE_CXX_STANDARD} standard.")

if (APPLE) 
    set(CMAKE_MACOSX_RPATH 1)
endif()

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.10.0')
    message(FATAL_ERROR "pybind11 version must be >= '2.10.0'")
endif()
add_definitions(-DPYBIND11_VERSION="${pybind11_VERSION}")

find_package(GSL REQUIRED)
option(USE_WEFFCPP "Use -Weffc++ during compilation" OFF)
option(ENABLE_PROFILING "Compile to enable code profiling" OFF)
option(BUILD_PYTHON_UNIT_TESTS "Build C++ modules for unit tests" ON)
option(BUILD_CPP_UNIT_TESTS "Build C++ unit test suite" OFF)
option(BUILD_CPP_BENCHMARK "Build C++ benchmark program" OFF)
option(DISABLE_LTO "Disable link-time optimization (LTO)" OFF)
option(CMAKE_BUILD "Set to ON when building with cmake. Set to off from, e.g., setup.py." ON)
include_directories(BEFORE
    ${fwdpy11_SOURCE_DIR}/fwdpy11/headers
    ${fwdpy11_SOURCE_DIR}/fwdpy11/headers/fwdpp
    ${fwdpy11_SOURCE_DIR}/lib)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})

if (USE_WEFFCPP)
    add_compile_options(-Weffc++)
endif()

if (ENABLE_PROFILING)
    add_definitions(-DPYBIND11_NAMESPACE=pybind11)
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wextra -Weffc++ -Woverloaded-virtual -Wold-style-cast -Werror=effc++ -Werror=old-style-cast -Werror=overloaded-virtual -Werror=unused-parameter")
endif()

if (CMAKE_BUILD)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/fwdpy11)
endif()
add_subdirectory(lib)
add_subdirectory(rust)
add_subdirectory(cpp)

if(BUILD_PYTHON_UNIT_TESTS)
    add_subdirectory(tests_with_cpp)
endif(BUILD_PYTHON_UNIT_TESTS)
if(BUILD_CPP_UNIT_TESTS)
    enable_testing()
    add_subdirectory(cpptests)
endif(BUILD_CPP_UNIT_TESTS)
if(BUILD_CPP_BENCHMARK)
    enable_testing()
    add_subdirectory(cpp_neutral_benchmark)
endif(BUILD_CPP_BENCHMARK)
