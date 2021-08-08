enable_language(Fortran)

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB srcs tecplot_merger_program/tecplot_merger.f90)
set(
        prog
        ${srcs}
)

# Example application target
add_executable(
        ${P_UTILITIES}
        "${prog}"
)

set_target_properties(
        ${P_UTILITIES}
        PROPERTIES
        OUTPUT_NAME "${PROJECT_NAME}"
)

target_link_libraries(
        ${P_UTILITIES}
        PRIVATE
        ${P_LIBRARY}
)


# Only install executable in main project
if(NOT is-subprojects)
    install(
            TARGETS
            ${P_UTILITIES}
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
