set(LIBRARY_SRC 
	${SRC_LIB}/src/boundary_conditions_class.f90
	${SRC_LIB}/src/boundary_type_class.f90
	${SRC_LIB}/src/chemical_properties_class.f90
	${SRC_LIB}/src/computational_domain_class.f90
	${SRC_LIB}/src/computational_mesh_class.f90
	${SRC_LIB}/src/data_io_class.f90
	${SRC_LIB}/src/data_manager_class.f90
	${SRC_LIB}/src/data_save_class.f90
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
)

add_library(${LIBRARY_NAME} STATIC ${LIBRARY_SRC})

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Make sure modules are exported
target_include_directories(${LIBRARY_NAME}
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)