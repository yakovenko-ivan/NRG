set(LIBRARY_SRC
    ${SRC_LIB}/src/boundary_conditions_class.f90
    ${SRC_LIB}/src/boundary_type_class.f90
    ${SRC_LIB}/src/chemical_properties_class.f90
    ${SRC_LIB}/src/computational_domain_class.f90
    ${SRC_LIB}/src/computational_mesh_class.f90
    ${SRC_LIB}/src/data_io_class.f90
    ${SRC_LIB}/src/data_manager_class.f90
    ${SRC_LIB}/src/field_pointers.f90
    ${SRC_LIB}/src/field_scalar_class.f90
    ${SRC_LIB}/src/field_tensor_class.f90
    ${SRC_LIB}/src/field_vector_class.f90
    ${SRC_LIB}/src/global_data.f90
    ${SRC_LIB}/src/kind_parameters.f90
    ${SRC_LIB}/src/mpi_communications_class.f90
    ${SRC_LIB}/src/post_processor_class.f90
    ${SRC_LIB}/src/post_processor_manager_class.f90
    ${SRC_LIB}/src/post_processor_operation_class.f90
    ${SRC_LIB}/src/solver_options_class.f90
    ${SRC_LIB}/src/supplementary_routines.f90
    ${SRC_LIB}/src/benchmarking_module.f90
    ${SRC_LIB}/src/thermophysical_properties_class.f90
    ${SRC_LIB}/src/run_control_class.f90

    # Format-independent full-field output layer.
    ${SRC_LIB}/src/field_output_types.f90
    ${SRC_LIB}/src/field_writer_class.f90
    ${SRC_LIB}/src/tecplot_writer_class.f90
    ${SRC_LIB}/src/cgns_writer_class.f90
    ${SRC_LIB}/src/data_save_class.f90
)

add_library(${LIBRARY_NAME} STATIC ${LIBRARY_SRC})

target_include_directories(${LIBRARY_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# -----------------------------------------------------------------------------
# Optional CGNS/HDF5 scientific field-output backend
# -----------------------------------------------------------------------------
option(NRG_ENABLE_CGNS "Enable CGNS/HDF5 full-field output" OFF)
option(NRG_USE_SYSTEM_CGNS
    "Use a preinstalled CGNS/HDF5 library instead of NRG's pinned bundled stack"
    OFF)

if(NRG_USE_SYSTEM_CGNS AND NOT NRG_ENABLE_CGNS)
    message(WARNING
        "NRG_USE_SYSTEM_CGNS=ON has no effect because NRG_ENABLE_CGNS=OFF")
endif()

if(NRG_ENABLE_CGNS)
    # Enable C at directory scope before configuring HDF5/CGNS. CMake does not
    # permit enable_language() from inside a function.
    enable_language(C)

    include(${CMAKE_SOURCE_DIR}/cmake/NRGCGNSDependencies.cmake)
    nrg_configure_cgns_dependency(_NRG_CGNS_TARGET)

    # Keep the C bridge in its own target so NRG's Fortran-only compiler flags
    # (/fpp, -ffree-line-length-512, etc.) never leak into C compilation.
    add_library(nrg_cgns_bridge STATIC
        ${SRC_LIB}/src/nrg_cgns_bridge.c
    )
    target_compile_features(nrg_cgns_bridge PRIVATE c_std_99)
    target_link_libraries(nrg_cgns_bridge PUBLIC ${_NRG_CGNS_TARGET})
    set_target_properties(nrg_cgns_bridge PROPERTIES
        FOLDER "Dependencies/CGNS-HDF5"
    )

    target_compile_definitions(${LIBRARY_NAME} PRIVATE NRG_CGNS)
    target_link_libraries(${LIBRARY_NAME} PUBLIC nrg_cgns_bridge)

    message(STATUS "CGNS/HDF5 output enabled via ${_NRG_CGNS_TARGET}")
else()
    message(STATUS "CGNS/HDF5 output disabled (NRG_ENABLE_CGNS=OFF)")
endif()
