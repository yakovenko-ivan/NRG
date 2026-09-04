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


# -----------------------------------------------------------------------------
# NRG source revision and package-library compatibility metadata
# -----------------------------------------------------------------------------
# Both package_interface and computing_module link package_library.  Generate a
# single build-info module here so setup generation and computation report the
# same provenance fields.  The committed package_library tree is the primary
# compatibility fingerprint; the dirty flag tells us whether local library
# changes make that committed tree insufficient to describe the actual build.
set(NRG_GIT_BRANCH "unknown")
set(NRG_GIT_HEAD "unknown")
set(NRG_PACKAGE_LIBRARY_TREE "unknown")
set(NRG_PACKAGE_LIBRARY_DIRTY "unknown")

find_package(Git QUIET)

if(GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                rev-parse --abbrev-ref HEAD
        RESULT_VARIABLE NRG_GIT_BRANCH_RESULT
        OUTPUT_VARIABLE NRG_GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT NRG_GIT_BRANCH_RESULT EQUAL 0)
        set(NRG_GIT_BRANCH "unknown")
    endif()

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                rev-parse HEAD
        RESULT_VARIABLE NRG_GIT_HEAD_RESULT
        OUTPUT_VARIABLE NRG_GIT_HEAD
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT NRG_GIT_HEAD_RESULT EQUAL 0)
        set(NRG_GIT_HEAD "unknown")
    endif()

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                rev-parse HEAD:package_library
        RESULT_VARIABLE NRG_PACKAGE_LIBRARY_TREE_RESULT
        OUTPUT_VARIABLE NRG_PACKAGE_LIBRARY_TREE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT NRG_PACKAGE_LIBRARY_TREE_RESULT EQUAL 0)
        set(NRG_PACKAGE_LIBRARY_TREE "unknown")
    endif()

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                status --porcelain --untracked-files=normal -- package_library
        RESULT_VARIABLE NRG_PACKAGE_LIBRARY_STATUS_RESULT
        OUTPUT_VARIABLE NRG_PACKAGE_LIBRARY_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NRG_PACKAGE_LIBRARY_STATUS_RESULT EQUAL 0)
        if(NRG_PACKAGE_LIBRARY_STATUS STREQUAL "")
            set(NRG_PACKAGE_LIBRARY_DIRTY "no")
        else()
            set(NRG_PACKAGE_LIBRARY_DIRTY "yes")
        endif()
    endif()

    # Reconfigure after branch switches and new commits so the embedded Git HEAD
    # follows the source revision used for a rebuilt executable.
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                rev-parse --git-path HEAD
        RESULT_VARIABLE NRG_GIT_HEAD_PATH_RESULT
        OUTPUT_VARIABLE NRG_GIT_HEAD_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NRG_GIT_HEAD_PATH_RESULT EQUAL 0)
        get_filename_component(
            NRG_GIT_HEAD_PATH "${NRG_GIT_HEAD_PATH}"
            ABSOLUTE BASE_DIR "${CMAKE_SOURCE_DIR}"
        )
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${NRG_GIT_HEAD_PATH}")
    endif()

    if(NOT NRG_GIT_BRANCH STREQUAL "HEAD" AND
       NOT NRG_GIT_BRANCH STREQUAL "unknown")
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                    rev-parse --git-path "refs/heads/${NRG_GIT_BRANCH}"
            RESULT_VARIABLE NRG_GIT_REF_PATH_RESULT
            OUTPUT_VARIABLE NRG_GIT_REF_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(NRG_GIT_REF_PATH_RESULT EQUAL 0)
            get_filename_component(
                NRG_GIT_REF_PATH "${NRG_GIT_REF_PATH}"
                ABSOLUTE BASE_DIR "${CMAKE_SOURCE_DIR}"
            )
            if(EXISTS "${NRG_GIT_REF_PATH}")
                set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                    "${NRG_GIT_REF_PATH}")
            endif()
        endif()
    endif()

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_SOURCE_DIR}"
                rev-parse --git-path packed-refs
        RESULT_VARIABLE NRG_GIT_PACKED_REFS_RESULT
        OUTPUT_VARIABLE NRG_GIT_PACKED_REFS_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NRG_GIT_PACKED_REFS_RESULT EQUAL 0)
        get_filename_component(
            NRG_GIT_PACKED_REFS_PATH "${NRG_GIT_PACKED_REFS_PATH}"
            ABSOLUTE BASE_DIR "${CMAKE_SOURCE_DIR}"
        )
        if(EXISTS "${NRG_GIT_PACKED_REFS_PATH}")
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                "${NRG_GIT_PACKED_REFS_PATH}")
        endif()
    endif()
endif()

message(STATUS "NRG Git branch: ${NRG_GIT_BRANCH}")
message(STATUS "NRG Git HEAD: ${NRG_GIT_HEAD}")
message(STATUS "NRG package_library tree: ${NRG_PACKAGE_LIBRARY_TREE}")
message(STATUS "NRG package_library modified: ${NRG_PACKAGE_LIBRARY_DIRTY}")

# Reconfigure if a package-library source changes locally so the dirty flag is
# refreshed before the affected library is rebuilt.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${LIBRARY_SRC}
    "${CMAKE_CURRENT_SOURCE_DIR}/src/nrg_build_info.f90.in")

string(REPLACE "'" "''" NRG_GIT_BRANCH_F90 "${NRG_GIT_BRANCH}")
string(REPLACE "'" "''" NRG_GIT_HEAD_F90 "${NRG_GIT_HEAD}")
string(REPLACE "'" "''" NRG_PACKAGE_LIBRARY_TREE_F90
    "${NRG_PACKAGE_LIBRARY_TREE}")
string(REPLACE "'" "''" NRG_PACKAGE_LIBRARY_DIRTY_F90
    "${NRG_PACKAGE_LIBRARY_DIRTY}")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(NRG_BUILD_INFO_SOURCE
    "${CMAKE_CURRENT_BINARY_DIR}/generated/nrg_build_info.f90")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/nrg_build_info.f90.in"
    "${NRG_BUILD_INFO_SOURCE}"
    @ONLY
)
list(APPEND LIBRARY_SRC "${NRG_BUILD_INFO_SOURCE}")

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
