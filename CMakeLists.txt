cmake_minimum_required(VERSION 3.13)

project(NRG VERSION 1.0
    DESCRIPTION "Numerical Reactive Gasdynamics JIHT RAS package"
    LANGUAGES Fortran
)

# -----------------------------------------------------------------------------
# Build configuration
# -----------------------------------------------------------------------------
get_property(NRG_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(NRG_MULTI_CONFIG)
    message(STATUS "Generator '${CMAKE_GENERATOR}' is multi-config")
    message(STATUS "Select the configuration at build time with --config <Debug|Release>")
else()
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
            "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel" FORCE)
        message(STATUS "CMAKE_BUILD_TYPE was not specified; defaulting to Debug")
    endif()

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

message(STATUS "Operating system: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Fortran compiler version: ${CMAKE_Fortran_COMPILER_VERSION}")

# -----------------------------------------------------------------------------
# Portable compiler and platform options
# -----------------------------------------------------------------------------
option(NRG_ENABLE_OPENMP "Build NRG with OpenMP support" ON)

# Interface target carrying all common compile definitions, compiler options,
# link options, and optional OpenMP requirements. Linking package_library to
# this target propagates the same configuration to all NRG executables.
add_library(nrg_build_options INTERFACE)

if(WIN32)
    target_compile_definitions(nrg_build_options INTERFACE WIN)
endif()

if(NRG_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(nrg_build_options INTERFACE OpenMP::OpenMP_Fortran)
    target_compile_definitions(nrg_build_options INTERFACE OMP)
    message(STATUS "OpenMP enabled: ${OpenMP_Fortran_VERSION}")
else()
    message(STATUS "OpenMP disabled")
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # GNU Fortran uses the same option spelling on Linux and Windows/MinGW.
    target_compile_options(nrg_build_options INTERFACE
        -cpp
        -ffree-line-length-512
        "$<$<CONFIG:Debug>:-O0;-g;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>"
    )

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(IntelLLVM|Intel)$")
    # IntelLLVM is ifx. Intel is retained for compatibility with classic ifort.
    if(WIN32)
        target_compile_options(nrg_build_options INTERFACE
            /fpp
            /traceback
            "$<$<CONFIG:Debug>:/Od;/debug:full;/check:all;/fpe:0>"
            "$<$<CONFIG:Release>:/O2;/fpe:0>"
            "$<$<CONFIG:RelWithDebInfo>:/O2;/debug:full;/fpe:0>"
            "$<$<CONFIG:MinSizeRel>:/O1;/fpe:0>"
        )
        target_link_options(nrg_build_options INTERFACE
            /traceback
            "$<$<CONFIG:Debug>:/check:all;/fpe:0>"
            "$<$<NOT:$<CONFIG:Debug>>:/fpe:0>"
        )
    else()
        target_compile_options(nrg_build_options INTERFACE
            -fpp
            -traceback
            "$<$<CONFIG:Debug>:-O0;-g;-check;all;-fpe0>"
            "$<$<CONFIG:Release>:-O2;-fpe0>"
            "$<$<CONFIG:RelWithDebInfo>:-O2;-g;-fpe0>"
            "$<$<CONFIG:MinSizeRel>:-O1;-fpe0>"
        )
        target_link_options(nrg_build_options INTERFACE
            -traceback
            "$<$<CONFIG:Debug>:-check;all;-fpe0>"
            "$<$<NOT:$<CONFIG:Debug>>:-fpe0>"
        )
    endif()

else()
    message(FATAL_ERROR
        "Unsupported Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'. "
        "Supported compiler families are GNU (gfortran) and IntelLLVM (ifx).")
endif()

# -----------------------------------------------------------------------------
# Target names and source directories
# -----------------------------------------------------------------------------
set(INTERFACE package_interface)
set(COMPUTE computing_module)
set(UTILITY package_utilities)
set(LIBRARY_NAME package_library)

set(SRC ${CMAKE_SOURCE_DIR})
set(SRC_INTERFACE ${SRC}/package_interface)
set(SRC_UTILITY ${SRC}/package_utilities)
set(SRC_COMPUTE ${SRC}/computing_module)
set(SRC_LIB ${SRC}/package_library)

# Put all generated Fortran module files in one build-tree directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Keep build products out of the source tree.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# -----------------------------------------------------------------------------
# User-selectable package interface and utility sources
# -----------------------------------------------------------------------------
set(PACKAGE_INTERFACE_SOURCE
    "src/tests/classic_tests/1D_laminar_velocity.f90"
    CACHE STRING
    "Fortran source file relative to package_interface/ used for package_interface"
)
message(STATUS "PACKAGE_INTERFACE_SOURCE: ${PACKAGE_INTERFACE_SOURCE}")

set(PACKAGE_UTILITY_SOURCE
    "src/leading_point_velocity.f90"
    CACHE STRING
    "Fortran source file relative to package_utilities/ used for package_utilities"
)
message(STATUS "PACKAGE_UTILITY_SOURCE: ${PACKAGE_UTILITY_SOURCE}")

# -----------------------------------------------------------------------------
# NRG targets
# -----------------------------------------------------------------------------
add_subdirectory(${SRC_LIB})
add_subdirectory(${SRC_INTERFACE})
add_subdirectory(${SRC_UTILITY})
add_subdirectory(${SRC_COMPUTE})

# Propagate common compiler/platform/OpenMP requirements through the library to
# every executable that links package_library.
target_link_libraries(${LIBRARY_NAME} PUBLIC nrg_build_options)
