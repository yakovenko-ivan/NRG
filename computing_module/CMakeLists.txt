enable_language(Fortran)

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB srcs src/current_build/*.f90)
file(GLOB slatec src/third_party/slatec/*.f90)
set(
        prog
        ${srcs} ${slatec}
)

# Example application target
add_executable(${COMPUTING_MODULE} ${prog})

set_target_properties(
        ${COMPUTING_MODULE}
        PROPERTIES
        OUTPUT_NAME "${PROJECT_NAME}"
)

target_link_libraries(
        ${COMPUTING_MODULE}
        PRIVATE
        ${P_LIBRARY}
        ${MPI_Fortran_LIBRARIES}
)


# Only install executable in main project
if(NOT is-subprojects)
    install(
            TARGETS
            ${COMPUTING_MODULE}
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
