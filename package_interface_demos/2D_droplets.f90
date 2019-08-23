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
	type(field_scalar_cons)			,target		:: p, T, rho, rho_d01, rho_d02
	type(field_vector_cons)			,target		:: v, Y
	
	type(liquid_droplets_phase)					:: liquid_droplets
	
	type(post_processor_manager)				:: problem_post_proc_manager
	
	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	real(dkind)	,dimension(3)	:: cell_size
	integer		,dimension(3,2)	:: utter_loop 
	integer		,dimension(3,2)	:: observation_slice, summation_region
	integer						:: transducer_offset
	
	integer						:: log_unit	
	
	
	character(len=500)			:: initial_work_dir
	character(len=100)			:: work_dir
	character(len=10)			:: solver_name
	
	real(dkind)	:: domain_length, ignition_region
	real(dkind)	:: CFL_coeff
	real(dkind)	:: delta_x, offset
	real(dkind)	:: nu
	real(dkind)	:: high_pressure
	
	logical	:: stop_flag, chemical_flag

	integer	:: task1, task2, task3, task4
	integer	:: ierr
	
	ierr = getcwd(initial_work_dir)
	
	do task1 = 1, 1
	
		work_dir = ''
		
		select case(task1)
			case(1)
				work_dir = trim(work_dir) // 'cpm'
				solver_name = 'cpm'
		end select
		
		ierr = system('mkdir '// work_dir)

		ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
		ierr = chdir(work_dir)

		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
	
		delta_x			= 2.5e-05_dkind	
		offset			= 0.0_dkind
		
		domain_length	= 0.01_dkind
		ignition_region	= 20.0e-05_dkind
		
		
		problem_domain			= computational_domain_c(	dimensions 			=	1,						&
															cells_number 		=	(/1000,1,1/),			&
															coordinate_system	=	'cartesian'	,			&													
															lengths				=	reshape((/	0.0_dkind,0.0_dkind,0.0_dkind,					&
																								0.1_dkind,0.005_dkind,0.005_dkind/),(/3,2/)),	&
															axis_names			=	(/'x','y','z'/))
	
		problem_chemistry		= chemical_properties_c(chemical_mechanism_file_name	= 'KEROMNES.txt'	,	&
														default_enhanced_efficiencies	= 1.0_dkind			,	&
														E_act_units						= 'cal.mol')
					
		problem_thermophysics	= thermophysical_properties_c(	chemistry					= problem_chemistry			,	&
																thermo_data_file_name		= 'KEROMNES_THERMO.txt'		,	&
																transport_data_file_name	= 'KEROMNES_TRANSDATA.txt'	,	&
																molar_masses_data_file_name	= 'molar_masses.dat')
	
		problem_mpi_support		= mpi_communications_c(problem_domain)

		problem_data_manager	= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics)
					
			call problem_data_manager%create_boundary_conditions(	problem_boundaries				, &	
																	number_of_boundary_types	= 3	, &
																	default_boundary			= 1)
																
			call problem_data_manager%create_computational_mesh(problem_mesh)
			call problem_data_manager%create_scalar_field(p		,'pressure'				, 'p')
			call problem_data_manager%create_scalar_field(T		,'temperature'			, 'T')
			call problem_data_manager%create_scalar_field(rho	,'density'				, 'rho')
			call problem_data_manager%create_scalar_field(rho_d01	,'density_droplets01'	, 'rho_d01')	
			call problem_data_manager%create_scalar_field(rho_d02	,'density_droplets02'	, 'rho_d02')	
			
			call problem_data_manager%create_vector_field(v		,'velocity'						,'v'	,'spatial')
			call problem_data_manager%create_vector_field(Y		,'specie_molar_concentration'	,'Y'	,'chemical')	

		cell_size			= problem_mesh%get_cell_edges_length()
		utter_loop			= problem_domain%get_global_utter_cells_bounds()			
			
		transducer_offset	= 0.001 / cell_size(1)

		observation_slice		= utter_loop
		observation_slice(2,:)	= 1
		
		summation_region(:,1)	= (/-transducer_offset,1,1/)
		summation_region(:,2)	= (/transducer_offset,1,1/)
		
		problem_post_proc_manager = post_processor_manager_c(problem_data_manager,number_post_processors = 1)
		
			call problem_post_proc_manager%create_post_processor(problem_data_manager			,	&
																post_processor_name = "proc1"	,	&
																operations_number	= 8			,	&
																save_time			= 1.0_dkind	,	&
																save_time_units		= 'nanoseconds')
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'min_grad'		,operation_area = observation_slice	,grad_projection = 1)		
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'specie_production_chemistry(H2)'	,'sum'			,operation_area = summation_region	,operation_area_distance = (/0,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'specie_molar_concentration(H2)'	,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'specie_molar_concentration(H2)'	,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'density'							,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'density'							,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
				call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
												
				
			problem_data_save	= data_save_c(	problem_data_manager	,	&
											visible_fields_names	= (/'pressure'						,&
																		'temperature'					,&
																		'density'						,&
																		'velocity'						,&
																		'specie_molar_concentration'	,&
																		'velocity_of_sound'				,&
																		'temperature_particles01'		,&
																		'density_particles01'			,&
																		'volume_fraction_particles01'	,&
																		'velocity_particles01'	,&
																		'temperature_particles02'		,&
																		'density_particles02'			,&
																		'volume_fraction_particles02'	,&
																		'velocity_particles02'	,&																		
																		'full_energy'/) ,	&
											save_time				= 2.0_dkind		,	&
											save_time_units			= 'microseconds'		,	&
											save_format				= 'tecplot'			,	&
											data_save_folder		= 'data_save'		,	&
											debug_flag				= .false.)
	
		problem_data_io		= data_io_c(problem_data_manager				,	&
										check_time			= 5000.0_dkind	,	&
										check_time_units	= 'microseconds'	,	&
										output_time			= 14400.0_dkind	,	&
										data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************

		p%cells(:,:,:)			= 1.0_dkind*101325.0_dkind
		T%cells(:,:,:)			= 300.0_dkind
		Y%pr(2)%cells(:,:,:)	= 1.0_dkind	

		T%cells(500:520,:,:)	= 1500.0_dkind
		
		Y%pr(1)%cells(501:,:,:)			= 1.0_dkind					!# Hydrogen
		Y%pr(2)%cells(501:,:,:)			= 0.5_dkind					!# Oxygen
	!	Y%pr(3)%cells(1:400,1:100,:)	= 0.5_dkind* 3.762_dkind	!# Nitrogen

		
		Y%pr(2)%cells(1:500,:,:)	= 1.0_dkind				!# Oxygen
		Y%pr(3)%cells(1:500,:,:)	= 3.762_dkind			!# Nitrogen
		
		rho_d01%cells(:,:,:)		= 1e-08_dkind
		rho_d01%cells(530:,:,:)		= 100.0_dkind

		!***********************************************************************************************
	
		call problem_thermophysics%change_field_units_mole_to_dimless(Y)

		!****************************** Setting boundary conditions ************************************
		call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
														slip					= .false.	,	&
														conductive				= .false.	,	&
														wall_temperature		= 0.0_dkind	,	&
														wall_conductivity_ratio	= 0.0_dkind	,	&
														farfield_pressure		= 0.0_dkind	,	&
														farfield_temperature	= 0.0_dkind ,	&
														priority				= 1)	
														
		call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
														slip					= .true.	,	&
														conductive				= .false.	,	&
														wall_temperature		= 0.0_dkind	,	&
														wall_conductivity_ratio	= 0.0_dkind	,	&
														farfield_pressure		= 0.0_dkind	,	&
														farfield_temperature	= 0.0_dkind ,	&
														priority				= 2)															
										
		call problem_boundaries%create_boundary_type (	type_name				= 'outlet'	,	&
														farfield_pressure		= 101325.0_dkind	,	&
														farfield_temperature	= 300.0_dkind		,	&
														farfield_velocity		= 0.0_dkind			,	&	
														farfield_species_names	= (/'H2,O2,N2'/)	,	&	
														farfield_concentrations	= (/1.0_dkind, nu, nu * 3.762_dkind/)	,	&
														priority				= 3)
															
		problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 3
		problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 3 	

	!	problem_boundaries%bc_markers(:,utter_loop(2,1),:)	= 2
	!	problem_boundaries%bc_markers(:,utter_loop(2,2),:)	= 1 
		
	!	problem_boundaries%bc_markers(1:400,101:,:)					= 1
	!	problem_boundaries%bc_markers(401:600,utter_loop(2,2),:)		= 3
		
		
		!***********************************************************************************************
	
		problem_solver_options = solver_options_c(	solver_name					= solver_name		, &
													hydrodynamics_flag			= .true.			, &
													heat_transfer_flag			= .true.			, &
													molecular_diffusion_flag	= .true.			, &
													viscosity_flag				= .true.			, &
													chemical_reaction_flag		= .true.			, &
													additional_droplets_phases	= 2					, &
													CFL_flag					= .true.			, &
													CFL_coefficient				= 0.75_dkind		, &
													initial_time_step			= 1.0e-06_dkind)
		
		liquid_droplets%diameter				= 0.1e-03_dkind
		liquid_droplets%material_density		= 1000.0_dkind
		liquid_droplets%material_heat_capacity	= 900.0_dkind
		liquid_droplets%material_latent_heat	= 2.0e+06_dkind
		liquid_droplets%combustible				= .false.										
		
	!	solid_particles%radii					= 5.0e-06_dkind
	!	solid_particles%material_conductivity	= 1000.0_dkind
	!	solid_particles%material_density		= 1000.0_dkind
		
		call problem_solver_options%create_additional_phase("liquid droplets",liquid_droplets_parameters=liquid_droplets)
													
		liquid_droplets%diameter				= 0.2e-03_dkind
		liquid_droplets%material_density		= 2000.0_dkind
		liquid_droplets%material_heat_capacity	= 100.0_dkind
		liquid_droplets%material_latent_heat	= 5.0e+06_dkind
		liquid_droplets%combustible				= .true.
													
	!	solid_particles%radii					= 5.0e-06_dkind
	!	solid_particles%material_conductivity	= 1000.0_dkind
	!	solid_particles%material_density		= 1000.0_dkind
		
		call problem_solver_options%create_additional_phase("liquid droplets",liquid_droplets_parameters=liquid_droplets)		
		
		!****************************** Writing problem short description ******************************
	
		write(log_unit,'(A)') 'General description: '
		write(log_unit,'(A)') 'Main aim: testing.'
		write(log_unit,'(A)') 'Problem setup: Keromnes oxidation scheme.'
		write(log_unit,'(A)') 'Solver setup:'
		write(log_unit,'(A)') 'Validation and comparison: '
		write(log_unit,'(A)') '--------------------------------------------------------------------------'
	
		!***********************************************************************************************
														
		call problem_data_io%output_all_data(0.0_dkind,stop_flag,make_output = .true.)												
		call problem_data_save%save_all_data(0.0_dkind,stop_flag,make_save = .true.)
	
		ierr = chdir(initial_work_dir)
		
	end do
		
	continue
end program