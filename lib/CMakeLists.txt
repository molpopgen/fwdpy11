add_subdirectory(corrosion)
set(DISCRETE_DEMOGRAPHY_SOURCES discrete_demography/SetDemeSize.cc)

set(MUTATION_DOMINANCE_SOURCES mutation_dominance/dependency_injection.cc)

set(ALL_SOURCES ${DISCRETE_DEMOGRAPHY_SOURCES}
    ${MUTATION_DOMINANCE_SOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Set up the rust build deps
get_filename_component(DEMES_FORWARD_HEADER_LOCATION ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY CACHE)
corrosion_import_crate(MANIFEST_PATH rust/demes-forward-capi/Cargo.toml)
add_custom_target(header DEPENDS ${DEMES_FORWARD_HEADER_LOCATION}/lib/demes_forward.h)

add_library(fwdpy11core STATIC ${ALL_SOURCES})
add_dependencies(fwdpy11core cargo-build_demes-forward-capi header)
target_link_directories(fwdpy11core PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(fwdpy11core PUBLIC demes_forward_capi STATIC)

# Run cbindgen manually to generate the demes header
add_custom_command(OUTPUT ${DEMES_FORWARD_HEADER_LOCATION}/lib/demes_forward.h COMMAND cbindgen -l C --cpp-compat -o
    ${DEMES_FORWARD_HEADER_LOCATION}/demes_forward.h ${CMAKE_SOURCE_DIR}/lib/rust/demes-forward-capi)
