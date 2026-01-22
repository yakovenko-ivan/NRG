program package_interface

	use ifport	
	
	use kind_parameters
	use global_data
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

	implicit none
	
	type(computational_domain)					:: problem_domain
	type(data_manager)							:: problem_data_manager
	type(mpi_communications)					:: problem_mpi_support

	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics

	type(solver_options)						:: problem_solver_options
	
	type(computational_mesh)		,target		:: problem_mesh
	type(boundary_conditions)		,target		:: problem_boundaries
	type(field_scalar_cons)			,target		:: p, T, rho, rho_d01
    type(field_vector_cons)			,target		:: v, Y
	
	type(particles_phase)					    :: droplets
	
	type(post_processor_manager)				:: problem_post_proc_manager
	
	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	real(dp)	,dimension(3)	:: cell_size
	integer		,dimension(3,2)	:: utter_loop 
	integer		,dimension(3,2)	:: observation_slice, summation_region
	integer						:: transducer_offset
	
	integer						:: log_unit	
	
	
	character(len=500)			:: initial_work_dir
	character(len=100)			:: work_dir
	character(len=20)			:: solver_name
	
	real(dp)	:: domain_length, ignition_region
	real(dp)	:: CFL_coeff
	real(dp)	:: delta_x, offset
	real(dp)	:: nu_air, nu_water
	real(dp)	:: ambient_temperature, ambient_pressure, p_crit, T_crit, saturation_pressure
    real(dp)    :: a1, a2, a3, a4, a5, a6
    real(dp)	:: droplets_diameter
    real(dp)	:: humidity, tau, max_humidity
	
	logical	:: stop_flag, chemical_flag

	integer	:: task1, task2, task3, task4
	integer	:: ierr
	
	ierr = getcwd(initial_work_dir)
	
	do task1 = 1, 1
	
        droplets_diameter	= 100.0e-06_dp
        ambient_temperature	= 400.0_dp
        ambient_pressure    = 2.0_dp*101325.0_dp
        nu_air              = 1.0_dp
        humidity            = 0.3_dp
        
		work_dir = ''
		
		select case(task1)
			case(1)
				work_dir = trim(work_dir) // 'fds_low_mach'
				solver_name = 'fds_low_mach'
		end select
		
		ierr = system('mkdir '// work_dir)

		ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
		ierr = chdir(work_dir)

		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
				
		problem_domain			= computational_domain_c(	dimensions 			=	3,						&
															cells_number 		=	(/10,10,10/),			&
															coordinate_system	=	'cartesian'	,			&													
															lengths				=	reshape((/	0.0_dp,0.0_dp,0.0_dp,					&
																								1.0_dp,1.0_dp,1.0_dp/),(/3,2/)),	&
															axis_names			=	(/'x','y','z'/))
	
		problem_chemistry		= chemical_properties_c(chemical_mechanism_file_name	= 'KEROMNES.txt'	,	&
														default_enhanced_efficiencies	= 1.0_dp			,	&
														E_act_units						= 'cal.mol')
					
		problem_thermophysics	= thermophysical_properties_c(	chemistry					= problem_chemistry			,	&
																thermo_data_file_name		= 'KEROMNES_THERMO.txt'		,	&
																transport_data_file_name	= 'KEROMNES_TRANSDATA.txt'	,	&
																molar_masses_data_file_name	= 'molar_masses.dat')
	
		problem_solver_options = solver_options_c(	solver_name					= solver_name		, &
													hydrodynamics_flag			= .true.			, &
													heat_transfer_flag			= .true.			, &
													molecular_diffusion_flag	= .true.			, &
													viscosity_flag				= .true.			, &
													chemical_reaction_flag		= .false.			, &
                                                    grav_acc                    = (/0.0_dp, 0.0_dp, 0.0_dp/)    , &
													additional_particles_phases	= 1				    , &
													CFL_flag					= .false.			, &
													CFL_coefficient				= 0.75_dp		, &
													initial_time_step			= 1.0e-04_dp)
        
		!# Water droplets
		droplets%diameter						= droplets_diameter
		droplets%material						= "H2O"
		droplets%material_density				= 1000.0_dp			
		droplets%material_heat_capacity			= 4181.73215473605_dp			
		droplets%material_latent_heat			= 2454208.88831280_dp		
		droplets%material_boiling_temperature	= 373.15_dp		
		droplets%evaporating         			= .true.	
        droplets%heating             			= .true.	
        droplets%inertial            			= .false.	
      
		call problem_solver_options%create_additional_phase(particles_parameters=droplets) 
		
        
        
		problem_mpi_support		= mpi_communications_c(problem_domain)

		problem_data_manager	= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics,problem_solver_options)
					
			call problem_data_manager%create_boundary_conditions(	problem_boundaries				, &	
																	number_of_boundary_types	= 1	, &
																	default_boundary			= 1)
																
			call problem_data_manager%create_computational_mesh(problem_mesh)
			call problem_data_manager%create_scalar_field(p		,'pressure'				, 'p')
			call problem_data_manager%create_scalar_field(T		,'temperature'			, 'T')
			call problem_data_manager%create_scalar_field(rho	,'density'				, 'rho')
			
			call problem_data_manager%create_vector_field(v		,'velocity'						,'v'	,'spatial')
			call problem_data_manager%create_vector_field(Y		,'specie_molar_concentration'	,'Y'	,'chemical')	

		cell_size			= problem_mesh%get_cell_edges_length()
		utter_loop			= problem_domain%get_global_utter_cells_bounds()			
				
		problem_post_proc_manager = post_processor_manager_c(problem_data_manager,number_post_processors = 0)
				
		problem_data_save	= data_save_c(	problem_data_manager	,	&
										visible_fields_names	= (/'pressure'						,&
																	'temperature'					,&
																	'density'						,&
																	'velocity'						,&
																	'specie_molar_concentration'	,&
																	'energy_production_droplets01'		,&
																	'concentration_production_droplets01'/), &
										save_time				= 10.0_dp		,	&
										save_time_units			= 'milliseconds'		,	&
										save_format				= 'tecplot'			,	&
										data_save_folder		= 'data_save'		,	&
										debug_flag				= .false.)
	
		problem_data_io		= data_io_c(problem_data_manager				,	&
										check_time			= 5000.0_dp	,	&
										check_time_units	= 'microseconds'	,	&
										data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************

		p%cells(:,:,:)			= ambient_pressure 
		T%cells(:,:,:)			= ambient_temperature
		
		Y%pr(2)%cells(:,:,:)	= nu_air * (1.0_dp)    !# Oxygen
		Y%pr(3)%cells(:,:,:)	= nu_air * (3.762_dp)	!# Nitrogen
        
        
        !# W. Wagner, A. Pruss doi: 10.1063/1.555926  
        T_crit = 647.096_dp
        p_crit = 22.064e06_dp
                     
        a1 =  -7.85951783_dp 
        a2 =   1.84408259_dp
        a3 = -11.7866497_dp 
        a4 =  22.6807411_dp 
        a5 = -15.9618719_dp
        a6 =   1.80122502_dp
        
        tau = 1.0_dp - ambient_temperature/T_crit
        
        saturation_pressure     = p_crit * exp(T_crit / ambient_temperature * (a1 * tau + a2 * tau**(1.5_dp) + a3 * tau**(3.0_dp) + a4 * tau**(3.5_dp) + a5 * tau**(4.0_dp) + a6 * tau**(7.5_dp)))
        
        max_humidity = ambient_pressure / saturation_pressure
        
        nu_water = humidity * nu_air / (max_humidity - humidity)	
        
        Y%pr(7)%cells(:,:,:)	= nu_water !# Water vapor
		
		!***********************************************************************************************
	
		call problem_thermophysics%change_field_units_mole_to_dimless(Y)

		!****************************** Setting boundary conditions ************************************
		call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
														slip					= .true.	,	&
														conductive				= .false.	,	&
														wall_temperature		= 0.0_dp	,	&
														wall_conductivity_ratio	= 0.0_dp	,	&
														farfield_pressure		= 0.0_dp	,	&
														farfield_temperature	= 0.0_dp ,	&
														priority				= 1)	
															
		problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 1
		problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 1 	

		problem_boundaries%bc_markers(:,utter_loop(2,1),:)	= 1
		problem_boundaries%bc_markers(:,utter_loop(2,2),:)	= 1
        
		problem_boundaries%bc_markers(:,:,utter_loop(3,1))	= 1
		problem_boundaries%bc_markers(:,:,utter_loop(3,2))	= 1 
		
		!****************************** Writing problem short description ******************************
	
		write(log_unit,'(A)') 'General description: '
		write(log_unit,'(A)') 'Main aim: testing.'
		write(log_unit,'(A)') 'Problem setup: Keromnes oxidation scheme.'
		write(log_unit,'(A)') 'Solver setup:'
		write(log_unit,'(A)') 'Validation and comparison: '
		write(log_unit,'(A)') '--------------------------------------------------------------------------'
	
		!***********************************************************************************************
														
		call problem_data_io%output_all_data(0.0_dp,stop_flag,make_output = .true.)												
		call problem_data_save%save_all_data(0.0_dp,stop_flag,make_save = .true.)
	
		ierr = chdir(initial_work_dir)
		
	end do
		
	continue
end program