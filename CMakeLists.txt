cmake_minimum_required(VERSION 3.15)

# Use CMake's first-class MSVC-runtime abstraction. This policy must be set
# before project() enables a language that targets the MSVC ABI.
cmake_policy(SET CMP0091 NEW)

# Select static or dynamic Microsoft/Intel multithreaded runtime libraries on
# Windows. For Intel Fortran this maps to /MT(/MTd) or /MD(/MDd), respectively.
# The option is defined before project() because CMAKE_MSVC_RUNTIME_LIBRARY must
# be established before targets/languages are created.
set(NRG_STATIC_RUNTIME ON CACHE BOOL
    "Link Intel/MSVC runtime libraries statically on Windows")

if(NRG_STATIC_RUNTIME)
    set(CMAKE_MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

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
    message(STATUS
        "Select the configuration at build time with "
        "--config <Debug|Release|RelWithDebInfo|MinSizeRel>")
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

# Interface target carrying all common compile definitions/options and optional
# OpenMP requirements. Linking package_library to this target propagates the
# same configuration to all NRG executables.
add_library(nrg_build_options INTERFACE)

if(WIN32)
    target_compile_definitions(nrg_build_options INTERFACE WIN)
endif()

if(NRG_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)

    # CMake's FindOpenMP publishes -Qiopenmp as both a compile option and an
    # INTERFACE_LINK_OPTION for IntelLLVM. With Visual Studio generators,
    # MSBuild performs the final link through the MSVC-compatible linker path.
    # With ifx 2023.1 this produces a malformed Intel link warning even though
    # the option is ignored and the executable links successfully.
    #
    # NRG uses host CPU OpenMP only, so preserve the imported target's compile
    # options and runtime-library requirements while suppressing only this
    # problematic IntelLLVM link-driver flag.
    if(WIN32
       AND CMAKE_GENERATOR MATCHES "^Visual Studio"
       AND CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set_property(TARGET OpenMP::OpenMP_Fortran
            PROPERTY INTERFACE_LINK_OPTIONS "")
        message(STATUS
            "IntelLLVM/Visual Studio: suppressed redundant OpenMP link-driver flag")
    endif()

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
    #
    # Traceback, checking and floating-point exception controls are compiler
    # options. Do not forward them through target_link_options(): with Visual
    # Studio the final linker is link.exe, which does not understand these
    # Intel Fortran options.
    if(WIN32)
        target_compile_options(nrg_build_options INTERFACE
            /fpp
            /traceback

            # Full runtime checking for correctness debugging, but suppress the
            # extremely noisy informational message emitted whenever the
            # compiler creates an argument temporary.
            "$<$<CONFIG:Debug>:/Od;/debug:full;/check:all,noarg_temp_created;/fpe:0>"

            # Production configuration.
            "$<$<CONFIG:Release>:/O2;/fpe:0>"

            # Release-like optimized code with symbols and traceback. Keep /O2
            # after /debug:full so optimization remains enabled. This is the
            # preferred configuration for optimization-sensitive failures.
            "$<$<CONFIG:RelWithDebInfo>:/debug:full;/O2;/fpe:0>"

            "$<$<CONFIG:MinSizeRel>:/O1;/fpe:0>"
        )

        if(NRG_STATIC_RUNTIME)
            message(STATUS
                "Intel Windows runtime linkage: static (/MT for non-Debug, /MTd for Debug)")
            if(NRG_ENABLE_OPENMP)
                message(STATUS
                    "Intel OpenMP runtime remains dynamic on Windows (libiomp5md.dll)")
            endif()
        else()
            message(STATUS
                "Intel Windows runtime linkage: dynamic (/MD for non-Debug, /MDd for Debug)")
        endif()
    else()
        target_compile_options(nrg_build_options INTERFACE
            -fpp
            -traceback
            "$<$<CONFIG:Debug>:-O0;-g;-check;all;-fpe0>"
            "$<$<CONFIG:Release>:-O2;-fpe0>"
            "$<$<CONFIG:RelWithDebInfo>:-g;-O2;-fpe0>"
            "$<$<CONFIG:MinSizeRel>:-O1;-fpe0>"
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

# Keep build products out of the source tree and separate configurations.
#
# The $<CONFIG> generator expression works for both generator families:
#   * single-config: CMAKE_BUILD_TYPE (Debug, Release, ...)
#   * multi-config:  configuration selected with --config
#
# Using the expression also prevents multi-config generators from appending a
# second configuration directory automatically.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/bin/$<CONFIG>")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/lib/$<CONFIG>")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/lib/$<CONFIG>")

# CMake creates the configuration-specific subdirectories during generation
# and build. Create only their common parents at configure time; otherwise a
# literal directory named '$<CONFIG>' would be created.
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

message(STATUS
    "Runtime output directory: ${CMAKE_BINARY_DIR}/bin/<CONFIG>")

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
