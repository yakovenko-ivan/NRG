program package_interface

	use IFPORT
	
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

	character(len=5000)			:: initial_work_dir
	character(len=100)			:: work_dir
	character(len=20)			:: solver_name		
	
	real(dkind)			,dimension(3)	:: cell_size
    real(dkind)			,dimension(2,2) :: cell_sizes_boundaries
    integer				,dimension(3)	:: cells_number
    real(dkind)			,dimension(3)   :: lengths
    real(dkind)			,dimension(:,:)		,allocatable :: reference_coordinates
    integer				,dimension(:)		,allocatable :: reference_coordinates_number
    real(dkind)			,dimension(:,:)		,allocatable :: reference_cell_sizes
    character(len=20)	,dimension(:)		,allocatable :: mesh_types	! uniform, linear, exponential
    integer                     :: dimensions
	integer		,dimension(3,2)	:: utter_loop, observation_slice
	integer		,dimension(3,2)	:: upstream_transducer, downstream_transducer, summation_region
	integer						:: log_unit	
    
    real(dkind)		:: cell_width, time_step_value, save_time_const, channel_radius, channel_length, ignition_radius, ignition_center_x, ignition_center_y
    real(dkind)		:: pressure_initial, temperature_initial, temperature_ignition, nitrogen_perc
    real(dkind)		:: sum_region_x, sum_region_y
	
	logical	:: stop_flag

	real	:: r
	
	integer	:: task1, task2, task3
    
    integer:: i, j, dim
	
	ierr = getcwd(initial_work_dir)
	
	do task1 = 3, 3
    do task2 = 1, 1
    do task3 = 1, 8
	
		work_dir = ''
        dimensions = 2
        
        time_step_value = 0.25_dkind * 10.0e-6_dkind / (2.0_dkind * 1250.0_dkind)
        channel_radius	= 5.0e-03_dkind 
        channel_length	= 5.0e-03_dkind  
        
        ignition_center_x	= channel_radius/2.0_dkind
        ignition_center_y	= channel_length/2.0_dkind
        ignition_radius		= 1.0e-3_dkind
        
        lengths(1)		 = channel_radius
        lengths(2)		 = channel_length
        lengths(3)		 = 1.0_dkind
	
		select case(task1)
			case(1)
				work_dir = trim(work_dir) // 'fds_low_mach'
				solver_name = 'fds_low_mach'
			case(2)
				work_dir = trim(work_dir) // 'cpm'
				solver_name = 'cpm'
            case(3)
				work_dir = trim(work_dir) // '2D_MESH_GENERATION_TEST'
				solver_name = 'CABARET'
        end select
        
		ierr = system('mkdir '// work_dir)
        
        select case(task2)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) // "H2_O2"
				nitrogen_perc			= 0.0_dkind
                save_time_const			= 10.0_dkind
                pressure_initial        = 20.0e03_dkind
				temperature_initial		= 300.0_dkind
				temperature_ignition	= 1500.0_dkind
        end select
		
		ierr = system('mkdir '// work_dir)
        
        select case(task3)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) // "uniform_100"
                cell_width = 100.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 2
                reference_coordinates_number(2) = 2
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)
                reference_coordinates(2, 2) = lengths(2)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
                
            case(2)
				work_dir = trim(work_dir) // trim(fold_sep) // "uniform_20"
                cell_width = 20.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 2
                reference_coordinates_number(2) = 2
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)
                reference_coordinates(2, 2) = lengths(2)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind

			case(3)
				work_dir = trim(work_dir) // trim(fold_sep) // "exp_200_20"
                cell_width = 200.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				mesh_types(1) = 'exponential'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 2
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/10.0_dkind
                
            case(4)
				work_dir = trim(work_dir) // trim(fold_sep) // "lin_100_20"
                cell_width = 100.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				mesh_types(1) = 'linear'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 2
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/5.0_dkind
                
            case(5)
				work_dir = trim(work_dir) // trim(fold_sep) // "exp_50_10_50"
                cell_width = 50.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				mesh_types(1) = 'exponential'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 3
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)/2.0_dkind
				reference_coordinates(1, 3) = lengths(1)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(1, 3)	= cell_width/1.0_dkind
                
            case(6)
				work_dir = trim(work_dir) // trim(fold_sep) // "lin_50_10_50"
                cell_width = 50.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'uniform'
				mesh_types(1) = 'linear'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 0
				reference_coordinates_number(1) = 3
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)/2.0_dkind
				reference_coordinates(1, 3) = lengths(1)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(1, 3)	= cell_width/1.0_dkind
                
            case(7)
				work_dir = trim(work_dir) // trim(fold_sep) // "lin_50_10_50_cross"
                cell_width = 50.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'linear'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 3
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)/2.0_dkind
				reference_coordinates(1, 3) = lengths(1)
				reference_coordinates(2, 2) = lengths(2)/2.0_dkind
				reference_coordinates(2, 3) = lengths(2)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(1, 3)	= cell_width/1.0_dkind
				reference_cell_sizes(2, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(2, 3)	= cell_width/1.0_dkind
                
            case(8)
				work_dir = trim(work_dir) // trim(fold_sep) // "exp_50_10_50_cross"
                cell_width = 50.0e-6_dkind
				allocate(mesh_types(dimensions))
				mesh_types = 'exponential'
				allocate(reference_coordinates_number(dimensions))
				reference_coordinates_number = 3
				allocate(reference_coordinates(dimensions,maxval(reference_coordinates_number)))
				reference_coordinates = 0.0_dkind
				reference_coordinates(1, 2) = lengths(1)/2.0_dkind
				reference_coordinates(1, 3) = lengths(1)
				reference_coordinates(2, 2) = lengths(2)/2.0_dkind
				reference_coordinates(2, 3) = lengths(2)
				allocate(reference_cell_sizes(dimensions,maxval(reference_coordinates_number)))
				reference_cell_sizes		= cell_width/1.0_dkind
				reference_cell_sizes(1, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(1, 3)	= cell_width/1.0_dkind
				reference_cell_sizes(2, 2)	= cell_width/5.0_dkind
				reference_cell_sizes(2, 3)	= cell_width/1.0_dkind
                
        end select
		
		ierr = system('mkdir '// work_dir)
		
		ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
		ierr = chdir(work_dir)		
		
		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
        
		
        problem_domain			= computational_domain_c(	dimensions 						=	dimensions,                                                             &									
                                                                coordinate_system				=	'cartesian'	,														&													
                                                                lengths							=	reshape((/	0.0_dkind, 0.000_dkind, 0.0_dkind,                          &
																												lengths(1), lengths(2), lengths(3)/),(/3,2/)),				&	 
                                                                axis_names						=	(/'x','y','z'/),														&
																mesh_types						=	mesh_types,                                                             &
																reference_coordinates			=	reference_coordinates,													&
                                                                reference_coordinates_number	=	reference_coordinates_number,											&
                                                                reference_cell_sizes			=	reference_cell_sizes)
       
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
																	number_of_boundary_types	= 4	, &
																	default_boundary			= 1)
																
			call problem_data_manager%create_computational_mesh(problem_mesh)
			call problem_data_manager%create_scalar_field(p		,'pressure'						, 'p')
			call problem_data_manager%create_scalar_field(T		,'temperature'					, 'T')
			call problem_data_manager%create_scalar_field(rho	,'density'						, 'rho')
	
			call problem_data_manager%create_vector_field(v		,'velocity'						,'v'		,'spatial')
			call problem_data_manager%create_vector_field(Y		,'specie_molar_concentration'	,'Y'		,'chemical')

			utter_loop			= problem_domain%get_global_utter_cells_bounds()			
		
			problem_post_proc_manager = post_processor_manager_c(problem_data_manager,number_post_processors = 0)
		
			problem_data_save	= data_save_c(	problem_data_manager	,	&
												visible_fields_names	= (/'pressure'						,&
																			!'pressure_dynamic'				,&
																			'temperature'					,&
																			'density'						,&
																			'velocity'						,&
																			'specie_molar_concentration'	,&
																			'velocity_of_sound'				,&
                                                                            !'full_energy'					,&
                                                                            'energy_production_chemistry'	,&
																			'velocity_production_viscosity' & 
																			/) ,	&
												save_time				= save_time_const		,	&
												save_time_units			= 'microseconds'		,	&
												save_format				= 'tecplot'			,	&
												data_save_folder		= 'data_save'		,	&
												debug_flag				= .false.)
	
			problem_data_io		= data_io_c(problem_data_manager				,	&
											check_time			= 1000.0_dkind	,	&
											check_time_units	= 'microseconds'	,	&
											output_time			= 40.0_dkind	,	&
											data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************

			T%cells(:,:,:)      = temperature_initial       !# Ambient temperature
			p%cells(:,:,:)      = pressure_initial			!# Ambient pressure
		
            Y%pr(2)%cells(:,:,:)	= 1.0_dkind         				! Oxygen moles
			Y%pr(1)%cells(:,:,:)	= 2.0_dkind							! Hydrogen moles
			Y%pr(3)%cells(:,:,:)	= nitrogen_perc * (1.0_dkind + 2.0_dkind) / (1.0_dkind - nitrogen_perc)     ! Nitrogen moles
            
            
			cell_size = problem_mesh%get_cell_edges_length_loc(utter_loop(1,1),utter_loop(2,1),1)
			sum_region_x = - 3.0_dkind * cell_size(1)/2.0_dkind
			do i = utter_loop(1,1), utter_loop(1,2)
				cell_size		= problem_mesh%get_cell_edges_length_loc(i,utter_loop(2,1),1)
				sum_region_y	= - 3.0_dkind * cell_size(2)/2.0_dkind
				do j = utter_loop(2,1), utter_loop(2,2)
					cell_size = problem_mesh%get_cell_edges_length_loc(i,j,1)
					if (j == utter_loop(2,1)) then
						sum_region_x = sum_region_x + cell_size(1)
					end if
					sum_region_y = sum_region_y + cell_size(2)
                
					! Ignition sone
					if ((sum_region_x - ignition_center_x)**2 + (sum_region_y - ignition_center_y)**2 <= ignition_radius**2)   then
						T%cells(i,j,:)	= temperature_ignition  		!# Ignition source
					end if
				end do
			end do

		!***********************************************************************************************
	
			call problem_thermophysics%change_field_units_mole_to_dimless(Y)

		!****************************** Setting boundary conditions ************************************
			call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
															slip					= .false.	,	&		
															conductive				= .false.	,	&		!# .true.		- for isothermal wall	, .false.	- for adiabatic wall
															wall_temperature		= 0.0_dkind	,	&		!# 300.0_dkind	- for isothermal wall	, 0.0_dkind	- for adiabatic wall
															wall_conductivity_ratio	= 0.0_dkind	,	&		!# 1.0_dkind	- for isothermal wall	, 0.0_dkind	- for adiabatic wall
															priority				= 1)	
            
            call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
															slip					= .true.	,	&		
															conductive				= .false.	,	&		!# .true.		- for isothermal wall	, .false.	- for adiabatic wall
															wall_temperature		= 0.0_dkind	,	&		!# 300.0_dkind	- for isothermal wall	, 0.0_dkind	- for adiabatic wall
															wall_conductivity_ratio	= 0.0_dkind	,	&		!# 1.0_dkind	- for isothermal wall	, 0.0_dkind	- for adiabatic wall
															priority				= 2)
			
			call problem_boundaries%create_boundary_type (	type_name				= 'outlet'	,	&
                                                            farfield_pressure		= pressure_initial	,	&
                                                            farfield_temperature	= temperature_initial		,	&
                                                            farfield_velocity		= 0.0_dkind			,	&	
                                                            farfield_species_names	= (/'O2','H2','N2'/)					,	&	
                                                            farfield_concentrations	= (/1.0_dkind, 2.0_dkind, (nitrogen_perc * (1.0_dkind + 2.0_dkind) / (1.0_dkind - nitrogen_perc))/)	,	&
                                                            priority				= 3)

															
			call problem_boundaries%create_boundary_type (	type_name				= 'inlet'	,	&
															farfield_pressure		= pressure_initial	,	&
															farfield_temperature	= temperature_initial		,	&
															farfield_velocity		= 0.00_dkind			,	&	
															farfield_species_names	= (/'O2','H2','N2'/)					,	&	
															farfield_concentrations	= (/1.0_dkind, 2.0_dkind, (nitrogen_perc * (1.0_dkind + 2.0_dkind) / (1.0_dkind - nitrogen_perc))/)	,	&
															priority				= 4)
																													
														
            problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 1
            problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 1
            problem_boundaries%bc_markers(:,utter_loop(2,1),:)	= 1
            problem_boundaries%bc_markers(:,utter_loop(2,2),:)	= 3
            
            
			!***********************************************************************************************
	
		problem_solver_options = solver_options_c(	solver_name					= solver_name, &
													hydrodynamics_flag			= .true.		, &
													heat_transfer_flag			= .true.		, &
													molecular_diffusion_flag	= .true.		, &
													viscosity_flag				= .true.		, &
													chemical_reaction_flag		= .true.		, &
													CFL_flag					= .true.		, &
													CFL_coefficient				= 0.25_dkind		, &
													initial_time_step			= time_step_value)													
														
														
		!****************************** Writing problem short description ******************************
	
		write(log_unit,'(A)') 'General description: '
		write(log_unit,'(A)') 'Main aim: '
		write(log_unit,'(A)') '			 '
		write(log_unit,'(A)') 'Problem setup:'
		write(log_unit,'(A)') '				 '
		write(log_unit,'(A)') '				 '
		write(log_unit,'(A)') 'Solver setup: '
		write(log_unit,'(A)') 'Validation and comparison:'
		write(log_unit,'(A)') '--------------------------------------------------------------------------'

		!***********************************************************************************************
														
		call problem_data_io%output_all_data(0.0_dkind,stop_flag,make_output	= .true.)												
		call problem_data_save%save_all_data(0.0_dkind,stop_flag,make_save		= .true.)
		
		ierr = chdir(initial_work_dir)
        
        deallocate(mesh_types)
        deallocate(reference_coordinates_number)
        deallocate(reference_coordinates)
        deallocate(reference_cell_sizes)
		
		continue
    end do
    end do
    end do
end program