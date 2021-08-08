cmake_minimum_required(VERSION 3.15)
enable_language(Fortran)

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB srcs src/*.f90)

set(
        MODULE_OUTPUT_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/include"
)

# Create example library target
add_library(${P_LIBRARY} ${srcs})

set_target_properties(
        ${P_LIBRARY}
        PROPERTIES
        POSITION_INDEPENDENT_CODE TRUE
        OUTPUT_NAME "${PROJECT_NAME}"
        VERSION "${PROJECT_VERSION}"
        SOVERSION "${PROJECT_VERSION_MAJOR}"
        Fortran_MODULE_DIRECTORY "${MODULE_OUTPUT_DIR}"
)
target_link_libraries(
        ${P_LIBRARY}
        PRIVATE
        "${lib-deps}"
)
target_include_directories(
        ${P_LIBRARY}
        INTERFACE
        $<BUILD_INTERFACE:${MODULE_OUTPUT_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Export targets for other projects
install(
        TARGETS
        ${P_LIBRARY}
        EXPORT
        "${PROJECT_NAME}-targets"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
        DIRECTORY
        "${MODULE_OUTPUT_DIR}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)