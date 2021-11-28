enable_language(Fortran)

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB srcs *.f90)

foreach (examplesourcefile ${srcs})
    get_filename_component(example_name ${examplesourcefile} NAME_WE)
    add_executable(${P_EXAMPLES}-${example_name} ${examplesourcefile})

    set_target_properties(
            ${P_EXAMPLES}-${example_name}
            PROPERTIES OUTPUT_NAME "${PROJECT_NAME}"
    )

    target_link_libraries(${P_EXAMPLES}-${example_name} PRIVATE ${P_LIBRARY})

    # Only install executable in main project
    if (NOT is-subprojects)
        install(
                TARGETS
                ${P_EXAMPLES}-${example_name}
                RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        )
    endif ()
endforeach (examplesourcefile ${srcs})