# Add the source files
# Use the PACKAGE_INTERFACE_SOURCE variable defined in the parent scope.
set(package_interface_src ${PACKAGE_INTERFACE_SOURCE})

# --- Extract base filename for output naming (e.g., "hydrogen_fill" from full path) ---
get_filename_component(PACKAGE_INTERFACE_BASENAME ${PACKAGE_INTERFACE_SOURCE} NAME_WE)

# --- Validate the source file path ---
# Construct the full expected path relative to this CMakeLists.txt
set(FULL_SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_INTERFACE_SOURCE}")
if(NOT EXISTS "${FULL_SOURCE_PATH}")
    # Provide a clear, fatal error message listing available options
    message(FATAL_ERROR "
    [NRG Configuration Error]
    The specified PACKAGE_INTERFACE_SOURCE was not found:
        '${PACKAGE_INTERFACE_SOURCE}'
    Expected at: ${FULL_SOURCE_PATH}

    Please ensure the file exists in the package_interface/ directory.
    Common examples are:
        - src/tests/classic_tests/1D_laminar_velocity.f90
        - src/tests/classic_tests/hydrogen_fill.f90
        - src/your_new_problem.f90

    Re-configure using:
        cmake .. -DPACKAGE_INTERFACE_SOURCE=<correct_path>
    ")
endif()
message(STATUS "Building package_interface from: ${PACKAGE_INTERFACE_SOURCE}")

# Define the executable in terms of the source files
add_executable(${INTERFACE} ${package_interface_src})

# Links the main executable with the created library
target_link_libraries(${INTERFACE} ${LIBRARY_NAME})

target_include_directories(${INTERFACE}
    PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# --- OPTIONAL: Give executable a descriptive name in output directory ---
# This will create e.g., "package_interface_hydrogen_fill" in build/bin/
set_target_properties(${INTERFACE} PROPERTIES OUTPUT_NAME "package_interface_${PACKAGE_INTERFACE_BASENAME}")