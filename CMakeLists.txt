cmake_minimum_required(VERSION 3.14)

project(
        "NRG"
        LANGUAGES "Fortran"
        VERSION "0.1"
)

# executable and libraries names
set(COMPUTING_MODULE "${PROJECT_NAME}-computing-module")
set(P_LIBRARY "${PROJECT_NAME}-package-library")
set(P_UTILITIES "${PROJECT_NAME}-package-utilities")
set(P_EXAMPLES "${PROJECT_NAME}-examples")

# Collect source of library module
add_subdirectory(package_library)
# Collect source of computing module
add_subdirectory(computing_module)
# Collect source of utilities module
add_subdirectory(package_utilities)
# Collect source of examples
add_subdirectory(package_interface_demos)


# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE ${P_LIBRARY})

install(
        TARGETS
        "${PROJECT_NAME}"
        EXPORT
        "${PROJECT_NAME}-targets"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)


# Install exported targets
install(
        EXPORT "${PROJECT_NAME}-targets"
        NAMESPACE
        "${PROJECT_NAME}::"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)


# select right compiler flugs
if (ENV MATCHES "DEBUG_OMP")
    message(STATUS "USING DEBUG_OMP COMPILER FLAGS")
    set(COMPUTING_FLAGS -g -traceback -fpp -DOMP -check all -fpe:0)
    set(P_LIBRARY_FLAGS -traceback -fpp -DOMP -check all -fpe:0 -warn all)
    set(P_UTILITIES_FLAGS -traceback -fpp -DOMP -check all -fpe:0 -warn all)
elseif (ENV MATCHES "DEBUG_MPI")
    message(STATUS "USING DEBUG_MPI COMPILER FLAGS")
    set(COMPUTING_FLAGS -g -traceback -fpp -Dmpi -check all -fpe:0)
    set(P_LIBRARY_FLAGS -traceback -fpp -Dmpi -check all -fpe:0 -warn all)
    set(P_UTILITIES_FLAGS -traceback -fpp -Dmpi -check all -fpe:0 -warn all)
elseif (ENV MATCHES "RELEASE_SEQUENTIAL")
    message(STATUS "USING RELEASE_SEQUENTIAL COMPILER FLAGS")
    set(COMPUTING_FLAGS -fpp -fpe:0)
    set(P_LIBRARY_FLAGS -fpp -fpe:0)
    set(P_UTILITIES_FLAGS -fpp -fpe:0)
elseif (ENV MATCHES "RELEASE_OMP")
    message(STATUS "USING RELEASE_OMP COMPILER FLAGS")
    set(COMPUTING_FLAGS -fpp -DOMP -qopenmp -fpe:0)
    set(P_LIBRARY_FLAGS -g -traceback -fpp -Dmpi -fpe:0)
    set(P_UTILITIES_FLAGS -g -traceback -fpp -Dmpi -fpe:0)
else () #DEBUG_SEQUENTIAL
    message(STATUS "USING DEBUG_SEQUENTIAL COMPILER FLAGS")
    set(COMPUTING_FLAGS -g -traceback -fpp -check all -fpe:0)
    set(P_LIBRARY_FLAGS -fpp -g -traceback)
    set(P_UTILITIES_FLAGS -g -traceback -fpp)
endif ()


if (NOT WIN32 AND ${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
    target_compile_options(${COMPUTING_MODULE} PRIVATE ${COMPUTING_FLAGS})
    target_compile_options(${P_LIBRARY} PRIVATE ${P_LIBRARY_FLAGS})
    target_compile_options(${P_UTILITIES} PRIVATE ${P_UTILITIES_FLAGS})
endif ()