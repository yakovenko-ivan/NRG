# Optional CGNS/HDF5 dependency setup for NRG.
#
# Public options are declared by package_library/CMakeLists.txt:
#   NRG_ENABLE_CGNS      - enable the CGNS/HDF5 field-output backend.
#   NRG_USE_SYSTEM_CGNS  - use an already installed CGNS package instead of
#                          the pinned bundled dependency stack.
#
# Bundled versions are intentionally pinned to combinations reviewed for NRG.
# CGNS 4.5.2 is built against HDF5 1.14.6 using only their C APIs.

function(nrg_configure_cgns_dependency out_target)
    if(NOT NRG_ENABLE_CGNS)
        set(${out_target} "" PARENT_SCOPE)
        return()
    endif()

    if(NRG_USE_SYSTEM_CGNS)
        _nrg_find_system_cgns(_nrg_cgns_target)
        set(${out_target} "${_nrg_cgns_target}" PARENT_SCOPE)
        return()
    endif()

    # CGNS 4.5.2 itself requires CMake 3.20.  Keep the older CMake requirement
    # available to dependency-free NRG builds, but diagnose bundled CGNS early.
    if(CMAKE_VERSION VERSION_LESS 3.20)
        message(FATAL_ERROR
            "Bundled CGNS/HDF5 output requires CMake 3.20 or newer. "
            "Update CMake, configure with -DNRG_USE_SYSTEM_CGNS=ON, or disable "
            "the backend with -DNRG_ENABLE_CGNS=OFF.")
    endif()

    include(FetchContent)

    # Visual Studio/Xcode solution folders are off by default on older CMake
    # releases.  Enable them explicitly so bundled third-party targets stay
    # grouped under Dependencies/CGNS-HDF5 rather than cluttering NRG's root.
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    set(NRG_BUNDLED_CGNS_VERSION "4.5.2")
    set(NRG_BUNDLED_HDF5_VERSION "1.14.6")

    message(STATUS
        "NRG CGNS backend: bundled dependency setup requested "
        "(HDF5 ${NRG_BUNDLED_HDF5_VERSION} + CGNS ${NRG_BUNDLED_CGNS_VERSION})")
    message(STATUS
        "NRG CGNS backend: first configuration downloads and configures the "
        "dependencies; this can take several minutes. Download progress follows.")

    # FetchContent is quiet by default on CMake versions using the legacy
    # population implementation.  Temporarily expose population output so the
    # first configure cannot look stalled.  Together with archive downloads,
    # CMake prints real '[download N% complete]' progress messages.
    _nrg_save_cache_value(FETCHCONTENT_QUIET _nrg_saved_fetchcontent_quiet)
    set(FETCHCONTENT_QUIET OFF CACHE BOOL
        "Show NRG bundled dependency download progress" FORCE)

    # Preserve the two generic cache switches that HDF5 uses so enabling CGNS
    # does not silently change unrelated targets in a parent/superbuild.
    _nrg_save_cache_value(BUILD_SHARED_LIBS _nrg_saved_build_shared)
    _nrg_save_cache_value(BUILD_TESTING _nrg_saved_build_testing)

    # ------------------------------------------------------------------ HDF5
    # CGNS only needs the HDF5 C API.  Compression, HL, language bindings,
    # examples, tools and parallel support are intentionally disabled.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
    set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Build tests" FORCE)
    set(ONLY_SHARED_LIBS OFF CACHE BOOL "Only build shared HDF5 libraries" FORCE)

    set(HDF5_EXTERNALLY_CONFIGURED ON CACHE BOOL "HDF5 embedded in NRG" FORCE)
    set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build HDF5 Fortran support" FORCE)
    set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "Build HDF5 C++ support" FORCE)
    set(HDF5_BUILD_JAVA OFF CACHE BOOL "Build HDF5 Java support" FORCE)
    set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "Build HDF5 high-level API" FORCE)
    set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 tools" FORCE)
    set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 examples" FORCE)
    set(HDF5_BUILD_GENERATORS OFF CACHE BOOL "Build HDF5 test generators" FORCE)
    set(HDF5_BUILD_PARALLEL_TOOLS OFF CACHE BOOL "Build HDF5 parallel tools" FORCE)
    set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "Enable parallel HDF5" FORCE)
    set(HDF5_ENABLE_THREADSAFE OFF CACHE BOOL "Enable thread-safe HDF5" FORCE)
    set(HDF5_ENABLE_Z_LIB_SUPPORT OFF CACHE BOOL "Enable HDF5 zlib filters" FORCE)
    set(HDF5_ENABLE_SZIP_SUPPORT OFF CACHE BOOL "Enable HDF5 SZIP filters" FORCE)
    set(HDF5_ENABLE_PLUGIN_SUPPORT OFF CACHE BOOL "Enable HDF5 filter plugins" FORCE)
    set(HDF5_NO_PACKAGES ON CACHE BOOL "Disable HDF5 packaging" FORCE)

    # Official HDF5 1.14.6 release archive.  Archive downloads are both faster
    # than a Git checkout and allow CMake to display percentage progress.
    # SHA256 is the published checksum for hdf5-1.14.6.tar.gz.
    message(STATUS
        "NRG CGNS backend [1/4]: fetching HDF5 ${NRG_BUNDLED_HDF5_VERSION} "
        "source archive (~39 MB)...")
    FetchContent_Declare(nrg_hdf5
        URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
        URL_HASH SHA256=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        DOWNLOAD_NO_PROGRESS FALSE
    )
    # Keep third-party targets grouped away from NRG's own projects in IDEs.
    # The directory is later marked EXCLUDE_FROM_ALL so only targets that are
    # actual dependencies of NRG remain relevant to the parent build.
    set(_nrg_saved_cmake_folder "${CMAKE_FOLDER}")
    set(CMAKE_FOLDER "Dependencies/CGNS-HDF5")
    FetchContent_MakeAvailable(nrg_hdf5)
    set(CMAKE_FOLDER "${_nrg_saved_cmake_folder}")

    if(DEFINED nrg_hdf5_SOURCE_DIR AND IS_DIRECTORY "${nrg_hdf5_SOURCE_DIR}")
        set_property(DIRECTORY "${nrg_hdf5_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
    endif()

    message(STATUS
        "NRG CGNS backend [2/4]: HDF5 ${NRG_BUNDLED_HDF5_VERSION} source ready and configured.")

    # CGNS 4.5.x calls find_package(HDF5) during its own configuration.  The
    # HDF5 1.14.x build-tree hdf5-config.cmake is generated with an install-like
    # relative target-file layout; when HDF5 is embedded through FetchContent,
    # that config can therefore look for hdf5-targets.cmake in the wrong parent
    # directory.  Do not consume that generated package.  Instead expose the
    # already-created hdf5-static target through a tiny compatibility package
    # containing exactly the HDF5 variables used by CGNS.
    if(NOT TARGET hdf5-static)
        message(FATAL_ERROR
            "Bundled HDF5 ${NRG_BUNDLED_HDF5_VERSION} configured, but target "
            "hdf5-static was not created.")
    endif()

    set(_nrg_hdf5_package_dir
        "${CMAKE_BINARY_DIR}/nrg-cmake-packages/hdf5")
    file(MAKE_DIRECTORY "${_nrg_hdf5_package_dir}")

    configure_file(
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/hdf5-buildtree/hdf5-config.cmake.in"
        "${_nrg_hdf5_package_dir}/hdf5-config.cmake"
        @ONLY
    )

    set(hdf5_DIR "${_nrg_hdf5_package_dir}" CACHE PATH
        "NRG compatibility package for bundled HDF5" FORCE)
    set(HDF5_DIR "${_nrg_hdf5_package_dir}" CACHE PATH
        "NRG compatibility package for bundled HDF5" FORCE)

    # ------------------------------------------------------------------ CGNS
    set(CGNS_ENABLE_HDF5 ON CACHE BOOL "Enable CGNS HDF5 backend" FORCE)
    set(CGNS_ENABLE_FORTRAN OFF CACHE BOOL "Build CGNS Fortran API" FORCE)
    set(CGNS_BUILD_SHARED OFF CACHE BOOL "Build shared CGNS library" FORCE)
    set(CGNS_USE_SHARED OFF CACHE BOOL "Link CGNS tools to shared library" FORCE)
    set(CGNS_ENABLE_TESTS OFF CACHE BOOL "Build CGNS tests" FORCE)
    set(CGNS_BUILD_TESTING OFF CACHE BOOL "Build CGNS CTest suite" FORCE)
    set(CGNS_BUILD_CGNSTOOLS OFF CACHE BOOL "Build CGNSTools GUI" FORCE)
    set(CGNS_ENABLE_PARALLEL OFF CACHE BOOL "Enable PCGNS" FORCE)
    set(HDF5_NEED_MPI OFF CACHE BOOL "CGNS HDF5 MPI dependency" FORCE)
    set(HDF5_NEED_ZLIB OFF CACHE BOOL "CGNS HDF5 zlib dependency" FORCE)
    set(HDF5_NEED_SZIP OFF CACHE BOOL "CGNS HDF5 SZIP dependency" FORCE)

    # CGNS 4.5.2 release archive.  The checksum corresponds to the upstream
    # GitHub tag archive and makes the dependency content reproducible.
    message(STATUS
        "NRG CGNS backend [3/4]: fetching CGNS ${NRG_BUNDLED_CGNS_VERSION} "
        "source archive (~2 MB)...")
    FetchContent_Declare(nrg_cgns
        URL https://github.com/CGNS/CGNS/archive/refs/tags/v4.5.2.tar.gz
        URL_HASH SHA256=95075e1fd0b51d97b1b96b73ebe03b1a551fbcc9cd2b2b6f487ccccedcff5964
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        DOWNLOAD_NO_PROGRESS FALSE
    )
    set(CMAKE_FOLDER "Dependencies/CGNS-HDF5")
    FetchContent_MakeAvailable(nrg_cgns)
    set(CMAKE_FOLDER "${_nrg_saved_cmake_folder}")

    if(DEFINED nrg_cgns_SOURCE_DIR AND IS_DIRECTORY "${nrg_cgns_SOURCE_DIR}")
        set_property(DIRECTORY "${nrg_cgns_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
    endif()

    # CGNS 4.5.2 adds its six command-line utilities unconditionally even when
    # CGNS_BUILD_CGNSTOOLS=OFF (that switch controls only the GUI package).
    # NRG does not use these utilities.  Exclude them from default builds and,
    # if a developer explicitly builds one, keep its executable in the
    # dependency tree rather than polluting NRG's bin/<CONFIG> directory.
    set(_nrg_cgns_cli_tools
        cgnsnames
        cgnscheck
        cgnsconvert
        cgnsdiff
        cgnslist
        cgnscompress
    )
    foreach(_nrg_tool IN LISTS _nrg_cgns_cli_tools)
        if(TARGET ${_nrg_tool})
            set_target_properties(${_nrg_tool} PROPERTIES
                EXCLUDE_FROM_ALL TRUE
                FOLDER "Dependencies/CGNS-HDF5/Optional tools"
                RUNTIME_OUTPUT_DIRECTORY "${nrg_cgns_BINARY_DIR}/tools-bin/$<CONFIG>"
            )
            if(CMAKE_GENERATOR MATCHES "Visual Studio")
                set_target_properties(${_nrg_tool} PROPERTIES
                    EXCLUDE_FROM_DEFAULT_BUILD TRUE
                )
            endif()
        endif()
    endforeach()

    message(STATUS
        "NRG CGNS backend [4/4]: CGNS ${NRG_BUNDLED_CGNS_VERSION} source ready and configured.")

    if(TARGET CGNS::cgns-static)
        set(_nrg_cgns_target CGNS::cgns-static)
    elseif(TARGET cgns_static)
        set(_nrg_cgns_target cgns_static)
    else()
        message(FATAL_ERROR
            "Bundled CGNS ${NRG_BUNDLED_CGNS_VERSION} configured, but its "
            "static CMake target was not created.")
    endif()

    # Restore generic parent-project switches.  HDF5/CGNS-specific cache values
    # remain cached deliberately, making the configured dependency stack visible
    # in `cmake -LA` and deterministic on subsequent reconfiguration.
    _nrg_restore_cache_value(BUILD_SHARED_LIBS "${_nrg_saved_build_shared}")
    _nrg_restore_cache_value(BUILD_TESTING "${_nrg_saved_build_testing}")
    _nrg_restore_cache_value(FETCHCONTENT_QUIET "${_nrg_saved_fetchcontent_quiet}")

    message(STATUS
        "NRG CGNS backend: bundled HDF5/CGNS dependency setup complete.")
    set(${out_target} "${_nrg_cgns_target}" PARENT_SCOPE)
endfunction()

function(_nrg_find_system_cgns out_target)
    set(_target "")

    find_package(cgns CONFIG QUIET)
    foreach(_candidate
            CGNS::cgns-static
            CGNS::cgns_static
            CGNS::cgns-shared
            CGNS::cgns_shared
            CGNS::CGNS)
        if(TARGET ${_candidate})
            set(_target ${_candidate})
            break()
        endif()
    endforeach()

    if(NOT _target)
        find_package(CGNS CONFIG QUIET)
        foreach(_candidate
                CGNS::cgns-static
                CGNS::cgns_static
                CGNS::cgns-shared
                CGNS::cgns_shared
                CGNS::CGNS)
            if(TARGET ${_candidate})
                set(_target ${_candidate})
                break()
            endif()
        endforeach()
    endif()

    if(NOT _target)
        # Fallback for installations without an exported CMake package.
        find_path(NRG_CGNS_INCLUDE_DIR NAMES cgnslib.h)
        find_library(NRG_CGNS_LIBRARY NAMES cgns cgnsdll)
        if(NOT NRG_CGNS_INCLUDE_DIR OR NOT NRG_CGNS_LIBRARY)
            message(FATAL_ERROR
                "NRG_USE_SYSTEM_CGNS=ON, but CGNS was not found. Set "
                "CMAKE_PREFIX_PATH/cgns_DIR, or provide NRG_CGNS_INCLUDE_DIR "
                "and NRG_CGNS_LIBRARY. The CGNS library must be built with HDF5.")
        endif()

        # A raw static CGNS library does not carry its HDF5 dependency metadata,
        # so recover that dependency explicitly for this compatibility fallback.
        find_package(HDF5 REQUIRED COMPONENTS C)

        if(NOT TARGET nrg_cgns_external)
            add_library(nrg_cgns_external UNKNOWN IMPORTED)
            set_target_properties(nrg_cgns_external PROPERTIES
                IMPORTED_LOCATION "${NRG_CGNS_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${NRG_CGNS_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${CMAKE_DL_LIBS}"
            )
        endif()
        set(_target nrg_cgns_external)
    endif()

    message(STATUS "NRG CGNS backend: using system target ${_target}")
    set(${out_target} "${_target}" PARENT_SCOPE)
endfunction()

# Store a cache variable as '<defined>|<type>|<value>'.
function(_nrg_save_cache_value variable out_value)
    if(DEFINED CACHE{${variable}})
        get_property(_type CACHE ${variable} PROPERTY TYPE)
        set(${out_value} "1|${_type}|${${variable}}" PARENT_SCOPE)
    else()
        set(${out_value} "0||" PARENT_SCOPE)
    endif()
endfunction()

function(_nrg_restore_cache_value variable saved)
    string(REPLACE "|" ";" _parts "${saved}")
    list(GET _parts 0 _was_defined)
    if(_was_defined)
        list(GET _parts 1 _type)
        list(GET _parts 2 _value)
        set(${variable} "${_value}" CACHE ${_type} "Restored by NRG dependency setup" FORCE)
    else()
        unset(${variable} CACHE)
    endif()
endfunction()
