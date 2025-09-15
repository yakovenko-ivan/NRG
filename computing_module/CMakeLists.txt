# Add the source files
set(computing_module_src 

	${SRC_COMPUTE}/src/third_party/slatec/d1mach.f90
	${SRC_COMPUTE}/src/third_party/slatec/daxpy.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddcor.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddcst.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddntl.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddntp.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddot.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddpsc.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddpst.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddriv3.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddscl.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddstp.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddzro.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgbfa.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgbsl.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgefa.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgesl.f90
	${SRC_COMPUTE}/src/third_party/slatec/dnrm2.f90
	${SRC_COMPUTE}/src/third_party/slatec/dscal.f90
	${SRC_COMPUTE}/src/third_party/slatec/fdump.f90
	${SRC_COMPUTE}/src/third_party/slatec/i1mach.f90
	${SRC_COMPUTE}/src/third_party/slatec/idamax.f90
	${SRC_COMPUTE}/src/third_party/slatec/j4save.f90
	${SRC_COMPUTE}/src/third_party/slatec/xercnt.f90
	${SRC_COMPUTE}/src/third_party/slatec/xerhlt.f90
	${SRC_COMPUTE}/src/third_party/slatec/xermsg.f90
	${SRC_COMPUTE}/src/third_party/slatec/xerprn.f90
	${SRC_COMPUTE}/src/third_party/slatec/xersve.f90
	${SRC_COMPUTE}/src/third_party/slatec/xgetua.f90

	${SRC_COMPUTE}/src/current_build/cabaret_solver.f90
	${SRC_COMPUTE}/src/current_build/chemical_kinetics_solver.f90
	${SRC_COMPUTE}/src/current_build/coarse_particles_method.f90
	${SRC_COMPUTE}/src/current_build/cpm_solver.f90
	${SRC_COMPUTE}/src/current_build/fds_low_mach_solver.f90
	${SRC_COMPUTE}/src/current_build/fickean_diffusion_solver.f90
	${SRC_COMPUTE}/src/current_build/fourier_heat_transfer_solver.f90
	${SRC_COMPUTE}/src/current_build/lagrangian_particles_solver.f90
	${SRC_COMPUTE}/src/current_build/main.f90
	${SRC_COMPUTE}/src/current_build/continuous_particles_solver.f90
	${SRC_COMPUTE}/src/current_build/riemann_solver.f90
	${SRC_COMPUTE}/src/current_build/table_approximated_real_gas.f90
	${SRC_COMPUTE}/src/current_build/viscosity_solver.f90
)

set(COMPUTE_BINARY_DIR ${SRC_COMPUTE}/bin)
file(MAKE_DIRECTORY ${COMPUTE_BINARY_DIR})


# Define the executable in terms of the source files
add_executable(${COMPUTE} ${computing_module_src})

# Group slatec files into a 'slatec' folder in Visual Studio
source_group("slatec" FILES 
	${SRC_COMPUTE}/src/third_party/slatec/d1mach.f90
	${SRC_COMPUTE}/src/third_party/slatec/daxpy.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddcor.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddcst.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddntl.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddntp.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddot.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddpsc.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddpst.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddriv3.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddscl.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddstp.f90
	${SRC_COMPUTE}/src/third_party/slatec/ddzro.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgbfa.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgbsl.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgefa.f90
	${SRC_COMPUTE}/src/third_party/slatec/dgesl.f90
	${SRC_COMPUTE}/src/third_party/slatec/dnrm2.f90
	${SRC_COMPUTE}/src/third_party/slatec/dscal.f90
	${SRC_COMPUTE}/src/third_party/slatec/fdump.f90
	${SRC_COMPUTE}/src/third_party/slatec/i1mach.f90
	${SRC_COMPUTE}/src/third_party/slatec/idamax.f90
	${SRC_COMPUTE}/src/third_party/slatec/j4save.f90
	${SRC_COMPUTE}/src/third_party/slatec/xercnt.f90
	${SRC_COMPUTE}/src/third_party/slatec/xerhlt.f90
	${SRC_COMPUTE}/src/third_party/slatec/xermsg.f90
	${SRC_COMPUTE}/src/third_party/slatec/xerprn.f90
	${SRC_COMPUTE}/src/third_party/slatec/xersve.f90
	${SRC_COMPUTE}/src/third_party/slatec/xgetua.f90
)


# Links the main executable with the created library
target_link_libraries(${COMPUTE} ${LIBRARY_NAME})

target_include_directories(${COMPUTE}
    PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

set_target_properties(${COMPUTE} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${COMPUTE_BINARY_DIR}/${CMAKE_BUILD_TYPE}
)
