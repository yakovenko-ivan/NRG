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

	type(computational_domain)					:: problem_domain
	type(data_manager)							:: problem_data_manager
	type(mpi_communications)					:: problem_mpi_support

	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics

	type(solver_options)						:: problem_solver_options
	
	type(computational_mesh)		,target		:: problem_mesh
	type(boundary_conditions)		,target		:: problem_boundaries
	type(field_scalar_cons)			,target		:: p, T, rho
	type(field_vector_cons)			,target		:: v, Y
	
	type(post_processor_manager)				:: problem_post_proc_manager
	
	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	real(dkind)	,dimension(3)	:: cell_size
	integer		,dimension(3,2)	:: utter_loop, observation_slice
	integer		,dimension(3,2)	:: upstream_transducer, downstream_transducer, summation_region
	integer						:: log_unit	
	
	character(len=500)			:: initial_work_dir
	character(len=100)			:: work_dir
	character(len=10)			:: solver_name

	real(dkind)	:: domain_length
	real(dkind)	:: CFL_coeff
	real(dkind)	:: delta_x
	real(dkind)	:: P_high
	
	logical	:: stop_flag

	integer	:: task1
	integer	:: ierr
	character*1  evalue
	ierr = getcwd(initial_work_dir)
	do task1 = 1, 3
	
		work_dir = ''
		
		select case(task1)
			case(1)
				work_dir = trim(work_dir) // 'cpm'
				solver_name = 'cpm'
			case(2)
				work_dir = trim(work_dir) // 'CABARET'
				solver_name = 'CABARET'
			case(3)
				work_dir = trim(work_dir) // 'CABARETM'
				solver_name = 'CABARETM'
		end select
!		character*1  evalue
		call getenv( 'USE_RIEMANN_MOD', evalue )
		print *, evalue
		work_dir =  trim(work_dir) // 'p1=20p2=1dx=4.0e-04CFL=0.5-'
		ierr = system('mkdir '// work_dir)

		delta_x	 = 4.0e-04_dkind
		CFL_coeff	 = 0.5_dkind
		ierr = system('cp -r .'//trim(fold_sep) // trim(task_setup_folder) // ' .' // trim(fold_sep) // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder))
		ierr = chdir('.' //trim(fold_sep) // work_dir)


		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')

		domain_length	= 0.5_dkind		
		
		problem_domain			= computational_domain_c(	dimensions 			=	1,						&
															cells_number 		=	(/int((domain_length)/delta_x),1,1/),			&
															coordinate_system	=	'cartesian'	,			&													
															lengths				=	reshape((/	0.0_dkind,0.0_dkind,0.0_dkind,					&
																								domain_length,0.005_dkind,0.005_dkind/),(/3,2/)),	&
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
																	number_of_boundary_types	= 1	, &
																	default_boundary			= 1)
																
			call problem_data_manager%create_computational_mesh(problem_mesh)
			call problem_data_manager%create_scalar_field(p		,'pressure'						, 'p')
			call problem_data_manager%create_scalar_field(T		,'temperature'					, 'T')
			call problem_data_manager%create_scalar_field(rho	,'density'						, 'rho')
	
			call problem_data_manager%create_vector_field(v		,'velocity'						,'v'		,'spatial')
			call problem_data_manager%create_vector_field(Y		,'specie_molar_concentration'	,'Y'		,'chemical')	
	
		utter_loop				= problem_domain%get_global_utter_cells_bounds()			
			
		problem_data_save	= data_save_c(	problem_data_manager	,	&
											visible_fields_names	= (/'pressure'						,&
																		'temperature'					,&
																		'density'						,&
																		'velocity'						,&
																		'specie_molar_concentration'	,&
																		'specie_production_chemistry'	,&
																		'energy_production_chemistry'	,&
																		'velocity_of_sound'				,&
																		'velocity_production_viscosity'	,&
																		'mixture_molar_concentration'	,&
																		'full_energy'/) ,	&
											save_time				= 5.0_dkind		,	&
											save_time_units			= 'microseconds'		,	&
											save_format				= 'tecplot'			,	&
											data_save_folder		= 'data_save'		,	&
											debug_flag				= .false.)
	
		problem_data_io		= data_io_c(problem_data_manager				,	&
										check_time			= 5000.0_dkind	,	&
										check_time_units	= 'nanoseconds'	,	&
										output_time			= 1400.0_dkind	,	&
										data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************
			! Ambient Temperature
			T%cells(:,:,:) = 300.0_dkind
			! Pressure in chamber
			p%cells(:,:,:) = 1.0_dkind*101325.0_dkind
			p%cells(:((0.25_dkind)/delta_x),:,:) = 20.0_dkind*101325.0_dkind
			! Velocity
!			v%pr(1)%cells(:,:,:) = 2000
!			v%pr(1)%cells(:((0.25_dkind)/delta_x),:,:) = 500
			! Velocity
!			v%pr(1)%cells(:,:,:) = 300.0_dkind
!			v%pr(1)%cells(:((0.25_dkind)/delta_x),:,:) = 2000.0_dkind
			! Air in entire tube
			Y%pr(1)%cells(:,:,:)	= 1.0_dkind					! Hydrogen
			!Y%pr(2)%cells(:,:,:)	= 0.5_dkind					! Oxygen
			!Y%pr(3)%cells(:,:,:)	= 1.881_dkind				! Nitrogen
			
			! Gaseous mixture in low-pressure chamber (Test mixture)
			!Y%pr(1)%cells(((0.5_dkind)/delta_x)+1:,:,:)	= 1.0_dkind					! Hydrogen
!			Y%pr(2)%cells(((0.5_dkind)/delta_x)+1:,:,:)	= 0.5_dkind					! Oxygen
			!Y%pr(3)%cells(((0.5_dkind)/delta_x)+1:,:,:)	= 1.881_dkind				! Nitrogen
			!
			! Gaseous mixture in high-pressure chamber (Driver gas)
			!Y%pr(4)%cells(:((0.5_dkind)/delta_x),:,:)	= 1.0_dkind					! Helium
			
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

		problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 1
		problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 1
	
		!***********************************************************************************************
	
		problem_solver_options = solver_options_c(	solver_name			= solver_name	, &
													hydrodynamics_flag			= .true.		, &
													heat_transfer_flag			= .true.		, &
													molecular_diffusion_flag	= .true.		, &
													viscosity_flag				= .true.		, &
													chemical_reaction_flag		= .false.		, &
													CFL_flag					= .true.		, &
													CFL_coefficient				= CFL_coeff		, &
													initial_time_step			= 1e-07_dkind)													
														
														
		!****************************** Writing problem short description ******************************
	
		write(log_unit,'(A)') 'General description: subsonic one-dimensional shock tube test with light driver (H2) and heavy driven (O2) gases.'
		write(log_unit,'(A)') 'Main aim: testing.'
		write(log_unit,'(A)') 'Problem setup: H2, O2, subsonic pressure drop across diaphragm 5atm - 1atm, 300K.'
		write(log_unit,'(A)') 'Solver setup: All processes including heat, diffusion and viscosity, varying timestep, initial time step 1e-07s.'
		write(log_unit,'(A)') 'Validation and comparison: analytical solution'
		write(log_unit,'(A)') '--------------------------------------------------------------------------'
	
		!***********************************************************************************************
														
		call problem_data_io%output_all_data(0.0_dkind,stop_flag,make_output = .true.)												
		call problem_data_save%save_all_data(0.0_dkind,stop_flag,make_save = .true.)
	
		ierr = chdir(initial_work_dir)

	end do	
		
	continue
end program