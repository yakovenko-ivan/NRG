# Add the source files
# Use the PACKAGE_UTILITY_SOURCE variable defined in the parent scope.
set(package_utilities_src ${PACKAGE_UTILITY_SOURCE})

# --- Extract base filename for output naming (e.g., "hydrogen_fill" from full path) ---
get_filename_component(PACKAGE_UTILITY_BASENAME ${PACKAGE_UTILITY_SOURCE} NAME_WE)

# --- Validate the source file path ---
# Construct the full expected path relative to this CMakeLists.txt
set(FULL_SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_UTILITY_SOURCE}")
if(NOT EXISTS "${FULL_SOURCE_PATH}")
    # Provide a clear, fatal error message listing available options
    message(FATAL_ERROR "
    [NRG Configuration Error]
    The specified PACKAGE_UTILITY_SOURCE was not found:
        '${PACKAGE_UTILITY_SOURCE}'
    Expected at: ${FULL_SOURCE_PATH}

    Please ensure the file exists in the package_interface/ directory.
    Common examples are:
        - src/leading_point_velocity.f90
        - src/output_merger.f90
        - src/your_new_problem.f90

    Re-configure using:
        cmake .. -DPACKAGE_UTILITY_SOURCE=<correct_path>
    ")
endif()
message(STATUS "Building package_utilities from: ${PACKAGE_UTILITY_SOURCE}")

# Define the executable in terms of the source files
add_executable(${UTILITY} ${package_utilities_src})

# Links the main executable with the created library
target_link_libraries(${UTILITY} ${LIBRARY_NAME})

target_include_directories(${UTILITY}
    PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# --- OPTIONAL: Give executable a descriptive name in output directory ---
# This will create e.g., "leading_point_velocity" in build/bin/
set_target_properties(${UTILITY} PROPERTIES OUTPUT_NAME "package_utilities_${PACKAGE_UTILITY_BASENAME}")