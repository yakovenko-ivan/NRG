# Add the source files
set(package_interface_src ${SRC_INTERFACE}/src/tests/classic_tests/1D_laminar_velocity.f90)

# Define the executable in terms of the source files
add_executable(${INTERFACE} ${package_interface_src})

# Links the main executable with the created library
target_link_libraries(${INTERFACE} ${LIBRARY_NAME})

target_include_directories(${INTERFACE}
    PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)
