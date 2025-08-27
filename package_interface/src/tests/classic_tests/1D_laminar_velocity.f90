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

	implicit none
	
	type(computational_domain)					:: problem_domain
	type(data_manager)							:: problem_data_manager
	type(mpi_communications)					:: problem_mpi_support

	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics

	type(solver_options)						:: problem_solver_options
	
	type(computational_mesh)		,target		:: problem_mesh
	type(boundary_conditions)		,target		:: problem_boundaries
	type(field_scalar_cons)			,target		:: p, T, rho, E_f, v_s
	type(field_vector_cons)			,target		:: v, Y
	
	type(post_processor_manager)				:: problem_post_proc_manager
	
	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	real(dp)	,dimension(3)	:: cell_size
	integer		,dimension(3,2)	:: utter_loop 
	integer		,dimension(3,2)	:: observation_slice, summation_region
	integer						:: transducer_offset
	integer						:: species_number
	
	integer						:: log_unit	
	
	integer						:: table_size, scaling_factor, n
	
	real(dp)	,allocatable	,dimension(:)		:: p_t,	T_t,	rho_t,	vx_t,	E_f_t,	gamma_t,	vs_t
	real(dp)	,allocatable	,dimension(:,:)		:: Y_t
	real(dp)	,allocatable	,dimension(:)		:: temp
	
	character(len=10)			:: string
	character(len=500)			:: initial_work_dir
	character(len=200)			:: work_dir, initials_file
	character(len=20)			:: solver_name, coordinate_system, setup
    character(len=30)			:: mech_file, thermo_file, transdata_file
	character(len=20)			:: str_r, str_i, str_e	
    
	real(dp)	:: domain_length, domain_width, ignition_region
	real(dp)	:: CFL_coeff
	real(dp)	:: delta_x, offset
	real(dp)	:: nu, X_H2, spec_summ
	
	logical	:: stop_flag

	
	integer	:: i, j, spec
	integer	:: task1, task2, task3, task4, task5
	integer	:: flamelet_pos
	integer	:: ierr, io_unit
	
	ierr = getcwd(initial_work_dir)
	
	do task1 = 1, 1
	do task2 = 1, 1
	do task3 = 1, 1
    do task4 = 9, 9
	do task5 = 2, 2
        
	
		work_dir = '1D_laminar_burning_velocity_test'
		
        ierr = system('mkdir '// work_dir)
	
		select case(task1)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) //'counter_flow'
				setup = 'counter_flow'
			case(2)
				work_dir = trim(work_dir) // trim(fold_sep) //'counter_flow_precInc'
				setup = 'counter_flow_precInc'                
			case(3)
				work_dir = trim(work_dir) // trim(fold_sep) //'near_wall'
				setup = 'near_wall'
		end select		
		
		ierr = system('mkdir '// work_dir)	
		
		select case(task2)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) // 'fds_low_mach'
				solver_name = 'fds_low_mach'
			case(2)
				work_dir = trim(work_dir) // trim(fold_sep) // 'cpm'
				solver_name = 'cpm'
			case(3)
				work_dir = trim(work_dir) // trim(fold_sep) // 'CABARET'
				solver_name = 'CABARET'				
        end select		
            
        ierr = system('mkdir '// work_dir)   
            
        select case(task3)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) // 'KEROMNES'
				mech_file		= 'KEROMNES.txt'
                thermo_file		= 'KEROMNES_THERMO.txt'
                transdata_file	= 'KEROMNES_TRANSDATA.txt'
			case(2)
				work_dir = trim(work_dir) // trim(fold_sep) // 'Tereza_FULL'
				mech_file		= 'Tereza_H2.txt'
                thermo_file		= 'Aramko20_Tereza_thermo.txt'
                transdata_file	= 'Aramko20_Tereza_transport.txt'
		end select
		
		ierr = system('mkdir '// work_dir)			

        X_H2		= task4 * 1.0_dp
        work_dir	= trim(work_dir) // trim(fold_sep) //  trim(str_r(X_H2)) //'_pcnt'
        nu			= (100.0 - X_H2) / X_H2 / 4.762_dp

		ierr = system('mkdir '// work_dir)			
		
		select case(task5)
			case(0)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_4.0e-04'
				delta_x	 = 4.0e-04_dp            
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_2.0e-04'
				delta_x	 = 2.0e-04_dp		
			case(2)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_1.0e-04'
				delta_x	 = 1.0e-04_dp
			case(3)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_5.0e-05'
				delta_x	 = 5.0e-05_dp
			case(4)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_2.5e-05'
				delta_x	 = 2.5e-05_dp
			case(5)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_1.25e-05'
				delta_x	 = 1.25e-05_dp
			case(6)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_6.25e-06'
				delta_x	 = 6.25e-06_dp		
		end select

		ierr = system('mkdir '// work_dir)			
		
		ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
		ierr = chdir(work_dir)		
		
		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
			
		domain_length	=	0.0256_dp
		
		problem_domain			= computational_domain_c(	dimensions 			=	1,						&
															cells_number 		=	(/int(domain_length/delta_x),1,1/),		&											
															coordinate_system	=	'cartesian'	,			&													
															lengths				=	reshape((/	0.0_dp,0.000_dp,0.0_dp,						&
																								domain_length,0.0025_dp,0.0025_dp/),(/3,2/)),	&	 
															axis_names			=	(/'x','y','z'/))
	
		problem_chemistry		= chemical_properties_c(chemical_mechanism_file_name	= mech_file			,	&
														default_enhanced_efficiencies	= 1.0_dp			,	&
														E_act_units						= 'cal.mol')
					
		problem_thermophysics	= thermophysical_properties_c(	chemistry					= problem_chemistry		,	&
																thermo_data_file_name		= thermo_file			,	&
																transport_data_file_name	= transdata_file		,	&
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
			
			cell_size			= problem_mesh%get_cell_edges_length()
			utter_loop			= problem_domain%get_global_utter_cells_bounds()			
			
			transducer_offset	= 0.005 / cell_size(1)
		
			observation_slice(:,1)		= (/1,1,1/)
			observation_slice(:,2)		= (/int(domain_length/delta_x),1,1/)			
			
			summation_region(:,1)	= (/-transducer_offset,1,1/)
			summation_region(:,2)	= (/transducer_offset,1,1/)		
		
			problem_post_proc_manager = post_processor_manager_c(problem_data_manager,number_post_processors = 0)

			!call problem_post_proc_manager%create_post_processor(problem_data_manager			,	&
			!													post_processor_name = "proc1"	,	&
			!													operations_number	= 7			,	&
			!													save_time			= 1.0_dp	,	&
			!													save_time_units		= 'microseconds')
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'min_grad'		,operation_area = observation_slice	,grad_projection = 1)		
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'pressure'						,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'pressure'						,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'density'							,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'density'							,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'transducer'	,operation_area_distance = (/transducer_offset,0,0/))
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,1,'temperature'						,'transducer'	,operation_area_distance = (/-transducer_offset,0,0/))
   !
			!call problem_post_proc_manager%create_post_processor(problem_data_manager			,	&
			!													post_processor_name = "proc2"	,	&
			!													operations_number	= 3			,	&
			!													save_time			= 1.0_dp	,	&
			!													save_time_units		= 'microseconds')
			!	call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,2,'pressure'						,'transducer'	,operation_area_distance = (/int(domain_lenght/delta_x),1,1/))
   !             call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,2,'temperature'						,'transducer'	,operation_area_distance = (/int(domain_lenght/delta_x),1,1/))
   !             call problem_post_proc_manager%create_post_processor_operation(problem_data_manager,2,'specie_molar_concentration(H2)'	,'transducer'	,operation_area_distance = (/int(domain_lenght/delta_x),1,1/))	               
		
		
			problem_data_save	= data_save_c(	problem_data_manager	,	&
												visible_fields_names	= (/'pressure'						,&
																			'pressure_dynamic'				,&
																			'temperature'					,&
																			'density'						,&
																			'velocity'						,&
																			'specie_molar_concentration'	,&
																			'velocity_of_sound'				,&
																			'velocity_production_viscosity'	,&
																			'mixture_cp'					,&
																			'thermal_conductivity'			,&
																			'viscosity'						,&
																			'energy_production_chemistry'	,&
																			'energy_production_diffusion'	&
																			/) ,	&
												save_time				= 100.0_dp		,	&
												save_time_units			= 'microseconds'	,	&
												save_format				= 'tecplot'			,	&
												data_save_folder		= 'data_save'		,	&
												debug_flag				= .false.)
	
			problem_data_io		= data_io_c(problem_data_manager					,	&
											check_time			= 5.0_dp			,	&
											check_time_units	= 'milliseconds'	,	&
											data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************

			T%cells(:,:,:)				= 300.0_dp					!# Ambient temperature
			p%cells(:,:,:)				= 1.0_dp*101325.0_dp	
		
			Y%pr(1)%cells(:,:,:)		= 1.0_dp
			Y%pr(2)%cells(:,:,:)		= nu
			Y%pr(3)%cells(:,:,:)		= nu * 3.762_dp
			
			select case(setup)
            case('counter_flow')
			!******** Ignition zone (counterflow) **********
				do i = utter_loop(1,1), utter_loop(1,2)
					if (i > int(3.0_dp*domain_length/delta_x/4.0_dp) - (0.005_dp/delta_x)) then
						T%cells(i,:,:)					= 2500.0_dp
                        if (i > int(3.0_dp*domain_length/delta_x/4.0_dp)) then
                            T%cells(i,:,:)				= 300.0_dp
							Y%pr(1)%cells(i,:,:)		= 0.0_dp				!# Hydrogen
							Y%pr(2)%cells(i,:,:)		= 0.0_dp				!# Oxygen
							Y%pr(3)%cells(i,:,:)		= nu * 3.762_dp		!# Nitrogen
							Y%pr(7)%cells(i,:,:)		= 1.0_dp	
                        end if
					end if			
                end do    
                
                call problem_thermophysics%change_field_units_mole_to_dimless(Y)
                
            case('counter_flow_precInc')
                
                initials_file	= trim(task_setup_folder) // trim(fold_sep) // 'table_initials' // trim(fold_sep) // trim('Tereza_FULL')// trim(fold_sep) // 'H2-Air_flamelet_' //  trim(str_r(X_H2)) //'_pcnt_dx_' // trim(str_e(2*delta_x)) // '.dat'
        
				open(newunit = io_unit, file = initials_file, status = 'old',form = 'formatted')								
				read(io_unit,*)	 string		
		
				table_size = 0
				do 
					read(io_unit,*,iostat = ierr)	string
					if(ierr /= 0) exit
					table_size = table_size + 1
				end do		
		
				table_size = table_size - 1
		
				species_number      = problem_chemistry%species_number
		
				allocate(Y_t(species_number,0:table_size-1), T_t(0:table_size-1),	vx_t(0:table_size-1))
				allocate(temp(species_number+4))
		
				rewind(io_unit)
				read(io_unit,*)	string
		
				do i = 0, table_size-1
					read(io_unit,*)	temp
					T_t(i)		=	temp(2)
					vx_t(i)		=	temp(4)
					do spec = 1, species_number
						Y_t(spec,i)	=	temp(spec+4)
					end do
				end do
        
				do i = utter_loop(1,1),utter_loop(1,2)

					flamelet_pos = int(2.0*domain_length/delta_x/4.0) 
                    
					if ((i >= flamelet_pos - (table_size - 1)).and.(i <= flamelet_pos + (table_size - 1))) then
                        if (mod(i - flamelet_pos,2) == 0) then 
							T%cells(i,1,1)					= T_t((i - flamelet_pos + (table_size-1))/2)
							v%pr(1)%cells(i,1,1)			= vx_t((i - flamelet_pos + (table_size-1))/2)
							do spec = 1, species_number
								Y%pr(spec)%cells(i,1,1)		= Y_t(spec,(i - flamelet_pos + (table_size-1))/2)
                            end do 
                        else
 							T%cells(i,1,1)					= 0.5 * (T_t((i - flamelet_pos + (table_size-1) + 1)/2) + T_t((i - flamelet_pos + (table_size-1) - 1)/2))
							v%pr(1)%cells(i,1,1)			= 0.5 * (vx_t((i - flamelet_pos + (table_size-1) + 1)/2) + vx_t((i - flamelet_pos + (table_size-1) - 1)/2))
							do spec = 1, species_number
								Y%pr(spec)%cells(i,1,1)		= 0.5 * (Y_t(spec,(i - flamelet_pos + (table_size-1) + 1)/2) + Y_t(spec,(i - flamelet_pos + (table_size-1) - 1)/2))
                            end do                            
                        
                        end if
					end if
				 !
					if(i<=flamelet_pos - (table_size - 1)) then
						T%cells(i,1,1)			= T_t(0)
						v%pr(1)%cells(i,1,1)	= vx_t(0)

                        spec_summ = 0.0
                        do spec = 1, species_number
							Y%pr(spec)%cells(i,1,1)	= Y_t(spec,0)
                            spec_summ = spec_summ + Y%pr(spec)%cells(i,1,1)
                        end do
              
                        do spec = 1, species_number 
							Y%pr(spec)%cells(i,1,1)	= Y%pr(spec)%cells(i,1,1) / spec_summ
                        end do    
					end if
		 
					if(i>flamelet_pos + (table_size - 1)) then
						T%cells(i,1,1)			= T_t(table_size-1)
						v%pr(1)%cells(i,1,1)	= vx_t(table_size-1)
						do spec = 1, species_number
							Y%pr(spec)%cells(i,1,1)	= Y_t(spec,table_size-1)
						end do
					end if
				end do
                
            case('near_wall')
			!******** Ignition zone (near the wall) **********
                j = 1
				do i = 1, floor(domain_length/delta_x)
                    if ((i < int(domain_length/delta_x/4))) then    
						T%cells(i,:,:)				= 2000.0_dp
						Y%pr(1)%cells(i,:,:)		= 0.0_dp				!# Hydrogen
						Y%pr(2)%cells(i,:,:)		= 0.0_dp				!# Oxygen
						Y%pr(3)%cells(i,:,:)		= nu * 3.762_dp		!# Nitrogen
						Y%pr(7)%cells(i,:,:)		= 1.0_dp
					end if
                end do		
                
                call problem_thermophysics%change_field_units_mole_to_dimless(Y)
            end select			
		!***********************************************************************************************



		!****************************** Setting boundary conditions ************************************
			call problem_boundaries%create_boundary_type (	type_name				= 'wall'	,	&
															slip					= .true.	,	&		
															conductive				= .false.	,	&		
															wall_temperature		= 0.0_dp	,	&		
															wall_conductivity_ratio	= 0.0_dp	,	&		
															priority				= 1)		

			select case(setup)
				case('counter_flow')															
					call problem_boundaries%create_boundary_type (	type_name				= 'outlet'							,	&
																	farfield_pressure		= 1.0_dp*101325.0_dp			,	&
																	farfield_temperature	= 1400.0_dp						,	&
																	farfield_velocity		= 0.0_dp							,	&	
																	farfield_species_names	= (/'H2O','N2'/)					,	&	
																	farfield_concentrations	= (/1.0_dp,nu * 3.762_dp/)	,	&
																	priority				= 2)
				case('near_wall')
					call problem_boundaries%create_boundary_type (	type_name				= 'outlet'							,	&
																	farfield_pressure		= 1.0_dp*101325.0_dp			,	&
																	farfield_temperature	= 300.0_dp						,	&
																	farfield_velocity		= 0.0_dp							,	&	
																	farfield_species_names	= (/'H2','O2','N2'/)				,	&	
																	farfield_concentrations	= (/1.0_dp,nu,nu * 3.762_dp/)	,	&
																	priority				= 2)	
               case('counter_flow_precInc')     
					call problem_boundaries%create_boundary_type (	type_name				= 'outlet'							,	&
																	farfield_pressure		= 1.0_dp*101325.0_dp			,	&
																	farfield_temperature	= T_t(table_size-1)					,	&
																	farfield_velocity		= vx_t(table_size-1)				,	&	
																	farfield_species_names	= (/'H2O','N2'/)					,	&	
																	farfield_concentrations	= (/1.0_dp,nu * 3.762_dp/)	,	&
																	priority				= 2)                    
            end select	
			
			select case(setup) 
                case('counter_flow')
					call problem_boundaries%create_boundary_type (	type_name				= 'inlet'							,	&
																	farfield_pressure		= 1.0_dp*101325.0_dp			,	&
																	farfield_temperature	= 300.0_dp						,	&
																	farfield_velocity		= 0.0_dp							,	&	
																	farfield_species_names	= (/'H2','O2','N2'/)				,	&	
																	farfield_concentrations	= (/1.0_dp,nu,nu * 3.762_dp/)	,	&
																	priority				= 3)
                case('counter_flow_precInc')
					call problem_boundaries%create_boundary_type (	type_name				= 'inlet'	,	&
																	farfield_pressure		= 1.0_dp*101325.0_dp	,	&
																	farfield_temperature	= T_t(0)				,	&
																	farfield_velocity		= vx_t(0)      ,	&	
																	farfield_species_names	= (/'H2','O2','N2'/)				,	&	
																	farfield_concentrations	= (/1.0_dp,nu,nu * 3.762_dp/)	,	&
																	priority				= 3)  
            end select
                    
															
			call problem_boundaries%create_boundary_type (	type_name				= 'symmetry_plane'					,	&
															priority				= 4)		

			select case(setup)
				case('counter_flow')
					problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 3	
					problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 2	
				case('counter_flow_precInc')
					problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 3	
					problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 2	                    
				case('near_wall')								
					problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 1
					problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 2	
			end select	
			
			!***********************************************************************************************
	
		problem_solver_options = solver_options_c(	solver_name					= solver_name, &
													hydrodynamics_flag			= .true.		, &
													heat_transfer_flag			= .true.		, &
													molecular_diffusion_flag	= .true.		, &
													viscosity_flag				= .true.		, &
													chemical_reaction_flag		= .true.		, &
													grav_acc					= (/0.0_dp, 0.0_dp, 0.0_dp/)	, &
													additional_droplets_phases	= 0				, &
													CFL_flag					= .true.		, &
													CFL_coefficient				= 0.25_dp	, &
													initial_time_step			= 1e-08_dp)													
										
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
														
		call problem_data_io%output_all_data(0.0_dp,stop_flag,make_output	= .true.)												
		call problem_data_save%save_all_data(0.0_dp,stop_flag,make_save		= .true.)
		
		ierr = chdir(initial_work_dir)	
		
		continue
		
	end do				
	end do		
	end do
	end do
	end do	

end program
    
character(len=20) function str_e(k)
	use kind_parameters

!   "Convert an real to string."
    real(dp), intent(in) :: k
    
    write (str_e, '(e10.3)') k    
    str_e = adjustl(str_e)
    
end function str_e
    
character(len=20) function str_r(k)
	use kind_parameters

!   "Convert an real to string."
    real(dp), intent(in) :: k
    
    write (str_r, '(f10.2)') k    
    str_r = adjustl(str_r)
    
end function str_r
    
character(len=20) function str_i(k)
	use kind_parameters

!   "Convert an integer to string."
    integer, intent(in) :: k
    
    write (str_i, '(I2)') k
    str_i = adjustl(str_i)
    
end function str_i
