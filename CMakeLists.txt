cmake_minimum_required(VERSION 3.10)

project(NRG VERSION 1.0
		DESCRIPTION "Numerical Reactive Gasdynmics JIHT RAS package"
		LANGUAGES Fortran)

enable_language(Fortran)

message(STATUS "Build Type is ${CMAKE_BUILD_TYPE}")

# Check whether the user has selected a BUILD_TYPE
IF ("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE "Debug")
	message(WARNING "Unable to deduce build type. Use -DCMAKE_BUILD_TYPE=<RELEASE/DEBUG>")
    message(STATUS "Build type not explicitly set. Trying to guess...")
	message(STATUS "Set to ${CMAKE_BUILD_TYPE}")
ELSE()
    message(STATUS "Build Type is ${CMAKE_BUILD_TYPE}")
ENDIF("${CMAKE_BUILD_TYPE}" STREQUAL "")

# if DEBUG requested, set the various debug options.
IF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG" OR ${CMAKE_BUILD_TYPE} MATCHES "Debug" OR ${CMAKE_GENERATOR} MATCHES "Visual Studio")
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
	message(NOTICE "Current Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ffree-line-length-512 -fcheck=all -fbacktrace -fbounds-check -cpp -DOMP -fopenmp")
    ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
		set(CMAKE_Fortran_FLAGS_DEBUG "/Od /fpe:0 /traceback /check:pointer /check:bounds /check:shape /check:uninit /check:format /check:output_conversion /check:stack /fpp /DWIN /DOMP /Qopenmp /libs:static /threads /c")
    ELSE()
        message(WARNING "Unable to determine Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    ENDIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message(STATUS "Debug flags for ${CMAKE_Fortran_COMPILER_ID}: ${CMAKE_Fortran_FLAGS_DEBUG}")
#ELSEIF
#	message(NOTICE "Current Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG" OR ${CMAKE_BUILD_TYPE} MATCHES "Debug" OR ${CMAKE_GENERATOR} MATCHES "Visual Studio")

IF(${CMAKE_BUILD_TYPE} MATCHES "RELEASE" OR ${CMAKE_BUILD_TYPE} MATCHES "Release" OR ${CMAKE_GENERATOR} MATCHES "Visual Studio")
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
	message(NOTICE "Current Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS_RELEASE "-ffree-line-length-512 -cpp -DOMP -fopenmp")
    ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
		set(CMAKE_Fortran_FLAGS_RELEASE "/nologo /debug:full /O2 /fpp /DWIN /DOMP /Qopenmp /fpe:0 /traceback /libs:static /threads /c")
    ELSE()
        message(WARNING "Unable to determine Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    ENDIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message(STATUS "Release flags for ${CMAKE_Fortran_COMPILER_ID}: ${CMAKE_Fortran_FLAGS_RELEASE}")
#ELSEIF
#	message(NOTICE "Current Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "RELEASE" OR ${CMAKE_BUILD_TYPE} MATCHES "Release" OR ${CMAKE_GENERATOR} MATCHES "Visual Studio")

# Define the executable name
set(INTERFACE package_interface)
# Define the executable name
set(COMPUTE computing_module)
# Name of source library
set(LIBRARY_NAME package_library)

# Source Directories
set(SRC ${CMAKE_SOURCE_DIR})
set(SRC_INTERFACE ${SRC}/package_interface)
set(SRC_COMPUTE ${SRC}/computing_module)
set(SRC_LIB ${SRC}/package_library)

# Set a global module output directory (optional but recommended)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# === NEW: Set global output directories for a clean build structure ===
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # For static libraries

# Ensure output directories exist
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
# === Define user-configurable source for package_interface ===
set(PACKAGE_INTERFACE_SOURCE "src/tests/classic_tests/1D_laminar_velocity.f90"
    CACHE STRING "Fortran source file (relative to package_interface/) to build the interface executable")
message(STATUS "PACKAGE_INTERFACE_SOURCE is set to: ${PACKAGE_INTERFACE_SOURCE}")

# Add subdirectories for the build process to look in
add_subdirectory(${SRC_LIB})
add_subdirectory(${SRC_INTERFACE})
add_subdirectory(${SRC_COMPUTE})