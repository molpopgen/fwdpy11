cmake_minimum_required(VERSION 2.8.12)
project(fwdpy11)

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.2.3')
    message(FATAL_ERROR "pybind11 version must be >= '2.2.3'")
endif()
add_subdirectory(fwdpy11)

