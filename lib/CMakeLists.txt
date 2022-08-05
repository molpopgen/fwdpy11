add_subdirectory(corrosion)
set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/SetDemeSize.cc)

set(MUTATION_DOMINANCE_SOURCES mutation_dominance/dependency_injection.cc)

set(ALL_SOURCES ${DISCRETE_DEMOGRAPHY_SOURCES}
    ${MUTATION_DOMINANCE_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Set up the rust build deps
get_filename_component(DEMES_FORWARD_HEADER_LOCATION ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY CACHE)
corrosion_import_crate(MANIFEST_PATH rust/demes-forward-capi/Cargo.toml)

add_library(fwdpy11core STATIC ${ALL_SOURCES})
