add_subdirectory(corrosion)
corrosion_import_crate(MANIFEST_PATH fp11rust/Cargo.toml)
add_custom_target(header DEPENDS ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h)
# NOTE: the build generates a warning
# that OpaqueForwardGraph is not available.
# This warning is a known false-positive.
# That symbol becomes visible at link-time and
# we provide a forward declaration in the relevant header file.
# We do not suppress warnings because we don't want to 
# suppress any that should be fixed.
add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h COMMAND cbindgen -l C --cpp-compat -o ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h ${CMAKE_SOURCE_DIR}/rust/fp11rust)
