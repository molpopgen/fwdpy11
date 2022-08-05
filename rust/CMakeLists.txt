add_subdirectory(corrosion)


# Set up the rust build deps
get_filename_component(DEMES_FORWARD_HEADER_LOCATION ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY CACHE)
corrosion_import_crate(MANIFEST_PATH demes-forward-capi/Cargo.toml)
add_custom_target(header DEPENDS ${DEMES_FORWARD_HEADER_LOCATION}/lib/demes_forward.h)
add_dependencies(cargo-build_demes-forward-capi header)

# Run cbindgen manually to generate the demes header
add_custom_command(OUTPUT ${DEMES_FORWARD_HEADER_LOCATION}/lib/demes_forward.h COMMAND cbindgen -l C --cpp-compat -o
    ${DEMES_FORWARD_HEADER_LOCATION}/lib/demes_forward.h ${CMAKE_SOURCE_DIR}/rust/demes-forward-capi)
