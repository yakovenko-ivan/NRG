program computing_module

	use IFPORT , only: SETENVQQ

	use global_data
	use kind_parameters
	use computational_domain_class
	use chemical_properties_class
	use thermophysical_properties_class
	use solver_options_class
	use computational_mesh_class
	use mpi_communications_class
	use data_manager_class
	use boundary_conditions_class
	use field_scalar_class
	use field_vector_class
	use data_save_class
	use data_io_class
	use post_processor_manager_class

	use cpm_solver_class	
	use cpm_low_mach_solver_class
	use cabaret_solver_class
	use fds_low_mach_solver_class
	
#ifdef mpi
	use MPI
#endif

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	type(computational_domain)					:: problem_domain
	type(mpi_communications)					:: problem_mpi_support
	type(computational_mesh)			,target	:: problem_mesh
	type(data_manager)							:: problem_manager
	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics
	type(cabaret_solver)						:: problem_cabaret_solver
	type(cpm_solver)							:: problem_cpm_solver
	type(cpm_low_mach_solver)					:: problem_cpm_low_mach_solver
	type(fds_solver)							:: problem_fds_solver

	type(solver_options)						:: problem_solver_options
	
	type(boundary_conditions)	,target			:: problem_boundaries
	type(field_scalar_cons)		,target			:: p, T, rho, e_i, E_f, mol_mix_conc
	type(field_vector_cons)		,target			:: v, Y

	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	type(post_processor_manager)				:: problem_post_proc_manager

	integer	:: c(8)
	integer	:: day, h, m , s

	integer		:: log_unit, mpi_io_unit

	integer		:: processor_rank

	integer		:: iter
	real(dkind)	:: calculation_time, time_step
	logical		:: stop_flag, precision_flag
	
	integer		:: error

#ifdef mpi
	call mpi_init(error)
#endif

#ifdef OMP
	call omp_set_num_threads(6)
#endif

	open(newunit = log_unit, file = problem_setup_log_file, status = 'old', form = 'formatted', position = 'append')

	problem_domain 			= computational_domain_c()

	problem_chemistry		= chemical_properties_c()
	problem_thermophysics	= thermophysical_properties_c(problem_chemistry)
	
	problem_mpi_support		= mpi_communications_c(problem_domain)

	problem_manager			= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics)

	problem_solver_options	= solver_options_c()
	
	call problem_manager%create_boundary_conditions(problem_boundaries)
	call problem_manager%create_computational_mesh(problem_mesh)
	call problem_manager%create_scalar_field(p				,'pressure'						,'P')
	call problem_manager%create_scalar_field(T 				,'temperature'					,'T')
	call problem_manager%create_scalar_field(rho			,'density'						,'rho')
	call problem_manager%create_scalar_field(e_i			,'internal_energy'				,'e_i')
	call problem_manager%create_scalar_field(E_f			,'full_energy'					,'E_f')
	call problem_manager%create_scalar_field(mol_mix_conc	,'mixture_molar_concentration'	,'mix_mol_conc')
	call problem_manager%create_vector_field(v				,'velocity'						,'v',	'spatial')
	call problem_manager%create_vector_field(Y				,'specie_molar_concentration '	,'Y',	'chemical')

!	problem_data_io				= data_io_c(problem_manager,calculation_time)

	select case(problem_solver_options%get_solver_name())
		case('cpm')
			problem_cpm_solver = cpm_solver_c(	problem_manager,  &
												problem_data_io			= problem_data_io		,	&
												problem_solver_options	= problem_solver_options)
		case('CABARET')
			problem_cabaret_solver = cabaret_solver_c(	problem_manager,  &
														problem_data_io			= problem_data_io		,	&
														problem_solver_options	= problem_solver_options)
		case('cpm_low_mach')											
			problem_cpm_low_mach_solver = cpm_low_mach_solver_c(	problem_manager,  &
																	problem_data_io			= problem_data_io		,	&
																	problem_solver_options	= problem_solver_options)
		case('fds_low_mach')											
			problem_fds_solver = fds_solver_c(	problem_manager,  &
												problem_data_io			= problem_data_io		,	&
												problem_solver_options	= problem_solver_options)																		
	end select

	problem_data_save			= data_save_c(problem_manager)
	problem_post_proc_manager	= post_processor_manager_c(problem_manager)	

	processor_rank = problem_domain%get_processor_rank()
	if (processor_rank == 0) then
		call problem_domain				%write_log(log_unit)													
		call problem_chemistry			%write_log(log_unit)		
		call problem_thermophysics		%write_log(log_unit)
		call problem_post_proc_manager	%write_log(log_unit)
		call problem_data_save			%write_log(log_unit)
		call problem_data_io			%write_log(log_unit)
		call problem_boundaries			%write_log(log_unit)
		call problem_solver_options		%write_log(log_unit)
	end if
	
	iter = 0
	stop_flag		= .false.
	precision_flag	= .false.
	do while(.not.stop_flag)
		iter = iter + 1

		select case(problem_solver_options%get_solver_name())
			case('cpm')
				call problem_cpm_solver%solve_problem()		
				calculation_time	= problem_cpm_solver%get_time()		
				time_step			= problem_cpm_solver%get_time_step()
			case('CABARET')
				call problem_cabaret_solver%solve_problem()
				calculation_time	= problem_cabaret_solver%get_time()
				time_step			= problem_cabaret_solver%get_time_step()
			case('cpm_low_mach')											
				call problem_cpm_low_mach_solver%solve_problem()		
				calculation_time	= problem_cpm_low_mach_solver%get_time()		
				time_step			= problem_cpm_low_mach_solver%get_time_step()		
			case('fds_low_mach')											
				call problem_fds_solver%solve_problem(iter,stop_flag)		
				calculation_time	= problem_fds_solver%get_time()		
				time_step			= problem_fds_solver%get_time_step()					
		end select
		
		if ((mod(iter-1,1) == 0).and.(processor_rank == 0)) then
			print *, ' Calculation time : ', calculation_time
			print *, ' Current time step : ', time_step
			print *, ' Amount of iterations : ', iter
  
			call date_and_time(values=c)
			day=c(3)
			h=c(5)
			m=c(6)
			s=c(7)
			print *, ' Current time = ', day,'  ',h,':',m,':',s
        end if
		
!		if ((precision_flag).and.(calculation_time > 165.0e-09_dkind)) then
!		if (time_step < 3.66e-08_dkind * problem_cpm_solver%get_CFL_coefficient()) then
!			call problem_data_save%set_save_time(1.0_dkind)
!			call problem_cpm_solver%set_CFL_coefficient(0.1_dkind)
!		else
!			call problem_cpm_solver%set_CFL_coefficient(0.75_dkind)
!		end if
	
!		if (calculation_time > 160.0e-04_dkind) then
!        if (calculation_time > 1.60e-04_dkind) then
!			call problem_data_save%set_save_time(1.0_dkind)
!        end if    
            
        call problem_data_io%output_all_data(calculation_time			,stop_flag)	
		call problem_post_proc_manager%process_data(calculation_time	,stop_flag)
		call problem_data_save%save_all_data(calculation_time			,stop_flag)

	end do

#ifdef mpi
	call mpi_finalize(error)
#endif

end program
