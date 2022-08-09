add_subdirectory(corrosion)
corrosion_import_crate(MANIFEST_PATH fp11rust/Cargo.toml)
add_custom_target(header DEPENDS ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h)
add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h COMMAND cbindgen -l C --cpp-compat -o ${CMAKE_SOURCE_DIR}/lib/core/internal/fp11rust_api.h ${CMAKE_SOURCE_DIR}/rust/fp11rust)
