cmake_minimum_required(VERSION 2.8.12)
project(fwdpy11)

set(CMAKE_CXX_EXTENSIONS OFF)
if(USECPP17)
    set(CMAKE_CXX_STANDARD 17) 
else()
    set(CMAKE_CXX_STANDARD 14) 
endif()
message(STATUS "Compiling with C++${CMAKE_CXX_STANDARD} standard.")

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.4.3')
    message(FATAL_ERROR "pybind11 version must be >= '2.4.3'")
endif()
add_definitions(-DPYBIND11_VERSION="${pybind11_VERSION}")

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

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wextra -Weffc++ -Woverloaded-virtual -Wold-style-cast -Werror=effc++ -Werror=old-style-cast -Werror=overloaded-virtual -Werror=unused-parameter")
endif()

add_subdirectory(fwdpy11)

if(BUILD_UNIT_TESTS)
    add_subdirectory(tests)
endif(BUILD_UNIT_TESTS)
