module cabaret_solver_class

	use kind_parameters
	use global_data
	use data_manager_class
	use data_io_class
	use computational_domain_class
	use computational_mesh_class
	use boundary_conditions_class
	use field_pointers
	use table_approximated_real_gas_class
	use thermophysical_properties_class	
	use chemical_properties_class

	use viscosity_solver_class
	use fickean_diffusion_solver_class
	use fourier_heat_transfer_solver_class
	use chemical_kinetics_solver_class

	use lagrangian_particles_solver_class
	use continuous_particles_solver_class
    
	use mpi_communications_class

	use benchmarking   
	use solver_options_class	
	
	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: cabaret_solver, cabaret_solver_c

	type(field_scalar_flow)	,target	::	rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new, T_f_new
	type(field_vector_flow)	,target	::	Y_f_new, v_f_new
    
    real(dp)	,dimension(:)	,allocatable	:: flame_front_coords
    integer	:: flame_loc_unit
    
#ifdef OMP
	integer(kind=omp_lock_kind)	,dimension(:,:,:)	,allocatable	:: lock
#endif	
	
    type(timer)     :: cabaret_timer
    type(timer)     :: cabaret_gas_dynamics_timer
    type(timer)     :: cabaret_eos_timer
    type(timer)     :: cabaret_chemistry_timer
    type(timer)     :: cabaret_diffusion_timer
    type(timer)     :: cabaret_heattransfer_timer
    type(timer)     :: cabaret_viscosity_timer

	type cabaret_solver
		logical			            :: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, hydrodynamics_flag, CFL_condition_flag
		real(dp)		            :: courant_fraction
		real(dp)		            :: time, time_step, initial_time_step
		real(dp)    , dimension(3)  :: g
		integer			:: additional_particles_phases_number
        
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(diffusion_solver)				:: diff_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(viscosity_solver)				:: viscosity_solver	
		type(table_approximated_real_gas)	:: state_eq

		type(lagrangian_particles_solver), dimension(:)	    ,allocatable	:: particles_solver			!# Lagrangian particles solver
		!type(continuous_particles_solver)   , dimension(:)	    ,allocatable	:: particles_solver			!# Continuum particles solver

        type(computational_domain)				:: domain
		type(mpi_communications)				:: mpi_support
		type(chemical_properties_pointer)		:: chem
		type(thermophysical_properties_pointer)	:: thermo
		type(computational_mesh_pointer)		:: mesh
		type(boundary_conditions_pointer)		:: boundary

		type(field_scalar_cons_pointer)	:: rho	, T	, p	, v_s, gamma, E_f	, e_i ,mol_mix_conc
		type(field_scalar_flow_pointer)	:: gamma_f_new, rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new, T_f_new
		
		type(field_scalar_cons_pointer)	:: E_f_prod_chem, E_f_prod_heat, E_f_prod_diff, E_f_prod_visc

		type(field_vector_cons_pointer)	:: v, Y	, v_prod_visc
		type(field_vector_flow_pointer)	:: v_f_new, Y_f_new
		
		type(field_vector_cons_pointer)	:: Y_prod_chem, Y_prod_diff

		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: rho_prod_particles, E_f_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	:: Y_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	:: v_prod_particles				!# Lagrangian particles solver        
        
		! Conservative variables
		real(dp) ,dimension(:,:,:)	,allocatable    :: rho_old, p_old, E_f_old, e_i_old, E_f_prod, rho_prod, v_s_old, gamma_old
		real(dp)	,dimension(:,:,:,:)	,allocatable	:: v_old, Y_old, v_prod, Y_prod

		! Flow variables
		real(dp) ,dimension(:,:,:,:,:)	,allocatable    :: v_f, Y_f
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: rho_f, p_f, e_i_f, E_f_f, v_s_f
		! Quasi invariants
		real(dp) ,dimension(:,:,:,:,:)		,allocatable    :: r_inv_corr, q_inv_corr, s_corr
        real(dp) ,dimension(:,:,:,:,:,:)		,allocatable    :: v_inv_corr, y_inv_corr
        
        
	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure	,private	:: apply_boundary_conditions_flow
		procedure				:: solve_problem
		procedure				:: solve_test_problem
		procedure				:: calculate_time_step
		procedure				:: get_time_step
		procedure				:: get_time
		procedure				:: set_CFL_coefficient
		procedure	,private	:: check_symmetry
        procedure	,private	:: if_stabilized
	end type

	interface	cabaret_solver_c
		module procedure	constructor
	end interface

contains

	type(cabaret_solver)	function constructor(manager,problem_data_io)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io

		real(dp)								:: calculation_time
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr		

		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr
		
        type(particles_phase)           :: particles_params
		integer							:: particles_phase_counter		
        
        character(len=40)       :: var_name
        
		real(dp)				:: spec_summ
		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds
		integer	,dimension(3,2)	:: flow_inner_loop, loop
		
		integer					:: dimensions, species_number
		integer					:: i, j, k, dim, dim1, spec
		real(dp)	,dimension(3)	:: cell_size

		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()

		species_number			= manager%chemistry%chem_ptr%species_number
		
		cell_size				= manager%computational_mesh_pointer%mesh_ptr%get_cell_edges_length()
		
		constructor%diffusion_flag		= manager%solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= manager%solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= manager%solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= manager%solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= manager%solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= manager%solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= manager%solver_options%get_CFL_condition_flag()
		constructor%sources_flag		= .false.
        
        constructor%g                       = manager%solver_options%get_grav_acc()

        constructor%additional_particles_phases_number	= manager%solver_options%get_additional_particles_phases_number()        
        
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%boundary%bc_ptr		=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
        
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr				=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr				=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'internal_energy')
		constructor%e_i%s_ptr				=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr		
		
		call manager%create_scalar_field(rho_f_new	,'density_flow'				,'rho_f_new')
		constructor%rho_f_new%s_ptr 	=> rho_f_new		
		call manager%create_scalar_field(p_f_new	,'pressure_flow'			,'p_f_new')
		constructor%p_f_new%s_ptr 		=> p_f_new
		call manager%create_scalar_field(e_i_f_new	,'internal_energy_flow'		,'e_i_f_new')
		constructor%e_i_f_new%s_ptr 	=> e_i_f_new
		call manager%create_scalar_field(E_f_f_new	,'full_energy_flow'			,'E_f_f_new')
		constructor%E_f_f_new%s_ptr 	=> E_f_f_new
		call manager%create_scalar_field(v_s_f_new	,'velocity_of_sound_flow'	,'v_s_f_new')
		constructor%v_s_f_new%s_ptr 	=> v_s_f_new
		call manager%create_scalar_field(T_f_new	,'temperature_flow'			,'T_f_new')
		constructor%T_f_new%s_ptr 		=> T_f_new
		
		call manager%create_vector_field(Y_f_new,'specie_molar_concentration_flow'	,'Y_f_new',	'chemical')
		constructor%Y_f_new%v_ptr => Y_f_new		
		call manager%create_vector_field(v_f_new,'velocity_flow'					,'v_f_new',	'spatial')
		constructor%v_f_new%v_ptr => v_f_new	
		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr	
		
		if (constructor%reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_ptr%v_ptr
		end if

		if (constructor%diffusion_flag) then
			constructor%diff_solver			= diffusion_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_diffusion')
			constructor%E_f_prod_diff%s_ptr			=> scal_ptr%s_ptr			
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_ptr%v_ptr
		end if

		if (constructor%heat_trans_flag) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr			=> scal_ptr%s_ptr
		end if

		if(constructor%viscosity_flag) then
			constructor%viscosity_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity')
			constructor%E_f_prod_visc%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_ptr%v_ptr
        end if			
		
        
		constructor%state_eq	=	table_approximated_real_gas_c(manager)        
        
		if(constructor%additional_particles_phases_number /= 0) then
			allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
			call constructor%particles_solver(1)%pre_constructor(constructor%additional_particles_phases_number)
			allocate(constructor%rho_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%Y_prod_particles(constructor%additional_particles_phases_number))
			do particles_phase_counter = 1, constructor%additional_particles_phases_number
				particles_params = manager%solver_options%get_particles_params(particles_phase_counter)
				constructor%particles_solver(particles_phase_counter)	= lagrangian_particles_solver_c(manager, particles_params, particles_phase_counter)		!# Lagrangian particles solver
!				constructor%particles_solver(particles_phase_counter)	= particles_solver_c(manager, particles_params, particles_phase_counter)                !# Continuum particles solver
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter						
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum particles solver								
				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum particles solver
				write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%Y_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr                
			end do		
		end if
	
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,'adiabatic_index_flow')
		constructor%gamma_f_new%s_ptr	=> scal_f_ptr%s_ptr				
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'adiabatic_index')
		constructor%gamma%s_ptr			=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_ptr%s_ptr		
		
		problem_data_io				= data_io_c(manager,calculation_time)
		
		if(problem_data_io%get_load_counter() /= 0) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%gamma)
			call problem_data_io%add_io_scalar_flow_field(constructor%rho_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%p_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%E_f_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%Y_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%v_f_new)
		end if

		call problem_data_io%input_all_data()

		if(problem_data_io%get_load_counter() == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%gamma)
			call problem_data_io%add_io_scalar_flow_field(constructor%rho_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%p_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%E_f_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%Y_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%v_f_new)
            
            if(constructor%additional_particles_phases_number /= 0) then
				do particles_phase_counter = 1, constructor%additional_particles_phases_number
					call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
				end do
			end if  
		end if		
					
		if(problem_data_io%get_load_counter() == 1) then
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
		else
			call constructor%state_eq%apply_state_equation()
			call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
		end if		
		
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%p%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%rho%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%E_f%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%T%s_ptr)
		call constructor%mpi_support%exchange_conservative_scalar_field(constructor%v_s%s_ptr)

		call constructor%mpi_support%exchange_conservative_vector_field(constructor%Y%v_ptr)
		call constructor%mpi_support%exchange_conservative_vector_field(constructor%v%v_ptr)

		call constructor%mpi_support%exchange_boundary_conditions_markers(constructor%boundary%bc_ptr)
		call constructor%mpi_support%exchange_mesh(constructor%mesh%mesh_ptr)

		allocate(constructor%rho_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
										
		allocate(constructor%p_old(		cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
										
		allocate(constructor%E_f_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
								
		allocate(constructor%e_i_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		allocate(constructor%E_f_prod(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		allocate(constructor%rho_prod(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		allocate(constructor%v_s_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		allocate(constructor%gamma_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		
		allocate(constructor%v_old(		dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%v_prod(	dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		allocate(constructor%Y_old(		species_number					, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
								
		allocate(constructor%Y_prod(	species_number					, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
										
			
										
		
		allocate(constructor%rho_f(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))
								
		allocate(constructor%p_f(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))
										
		allocate(constructor%E_f_f(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))	
										
		allocate(constructor%e_i_f(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))
										
		allocate(constructor%v_s_f(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))	
										
		allocate(constructor%v_f(		dimensions						, &
										dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))	
								
		allocate(constructor%Y_f(		species_number					, &
										dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))		

		allocate(constructor%r_inv_corr(			2						, &
										dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2))) 
        
		allocate(constructor%q_inv_corr(			2						, &
										dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))   
        
		allocate(constructor%s_corr(			2						, &
										dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))  
        
        
		allocate(constructor%v_inv_corr(			2						, &	
											dimensions						, &
											dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2))) 
        
		allocate(constructor%y_inv_corr(			2						, &	
											species_number					, &
											dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))   
        
        
        
		flow_inner_loop	= manager%domain%get_local_inner_faces_bounds()									
	
		
#ifdef OMP	
		allocate(lock(	flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
						flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
						flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		do k = flow_inner_loop(3,1),flow_inner_loop(3,2)
		do j = flow_inner_loop(2,1),flow_inner_loop(2,2)
		do i = flow_inner_loop(1,1),flow_inner_loop(1,2)		
			call omp_init_lock(lock(i,j,k))			
		end do
		end do
		end do
#endif											
										
		if (problem_data_io%get_load_counter() == 1) then

			constructor%p_f(:,:,:,:)	= 0.0_dp								
			constructor%rho_f(:,:,:,:)	= 0.0_dp
			constructor%Y_f(:,:,:,:,:)	= 0.0_dp
			constructor%v_f(:,:,:,:,:)	= 0.0_dp

			do dim = 1, dimensions		

				loop = flow_inner_loop

				do dim1 = 1,dimensions
					loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
				end do

				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2) 
				do i = loop(1,1),loop(1,2) 
						
			
					!	if ((constructor%boundary%bc_ptr%bc_markers(i,j,k) == 0).and.(constructor%boundary%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
							constructor%p_f(dim,i,j,k)		=	0.5_dp * (constructor%p%s_ptr%cells(i,j,k)	+ constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	
							constructor%rho_f(dim,i,j,k)	=	0.5_dp * (constructor%rho%s_ptr%cells(i,j,k)	+ constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	
							constructor%E_f_f(dim,i,j,k)	=	0.5_dp * (constructor%E_f%s_ptr%cells(i,j,k)	+ constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%v_s_f(dim,i,j,k)	=	0.5_dp * (constructor%v_s%s_ptr%cells(i,j,k)	+ constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							spec_summ = 0.0_dp
							do spec = 1, species_number
								constructor%Y_f(spec,dim,i,j,k) = 0.5_dp * (constructor%Y%v_ptr%pr(spec)%cells(i,j,k) + constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								spec_summ = spec_summ + max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp)
							end do

							do spec = 1,species_number
								constructor%Y_f(spec,dim,i,j,k) = max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp) / spec_summ
							end do

							do dim1 = 1, dimensions
								constructor%v_f(dim1,dim,i,j,k) = 0.5_dp * (constructor%v%v_ptr%pr(dim1)%cells(i,j,k) + constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) )
							end do

							if (constructor%boundary%bc_ptr%bc_markers(i,j,k) /= 0) then
								constructor%p_f(dim,i,j,k)		=	constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	
								constructor%rho_f(dim,i,j,k)	=	constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								constructor%E_f_f(dim,i,j,k)	=	constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								constructor%v_s_f(dim,i,j,k)	=	constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))

								spec_summ = 0.0_dp
								do spec = 1, species_number
									constructor%Y_f(spec,dim,i,j,k) =  constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
									spec_summ = spec_summ + max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp)
								end do

								do spec = 1,species_number
									constructor%Y_f(spec,dim,i,j,k) = max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp) / spec_summ
								end do

								do dim1 = 1, dimensions
									constructor%v_f(dim1,dim,i,j,k) =	constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end do
							end if

							if (constructor%boundary%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) /= 0) then
								constructor%p_f(dim,i,j,k)		=	constructor%p%s_ptr%cells(i,j,k)	
								constructor%rho_f(dim,i,j,k)	=	constructor%rho%s_ptr%cells(i,j,k)	
								constructor%E_f_f(dim,i,j,k)	=	constructor%E_f%s_ptr%cells(i,j,k)
								constructor%v_s_f(dim,i,j,k)	=	constructor%v_s%s_ptr%cells(i,j,k)
								
								spec_summ = 0.0_dp
								do spec = 1, species_number
									constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
									spec_summ = spec_summ + max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp)
								end do

								do spec = 1,species_number
									constructor%Y_f(spec,dim,i,j,k) = max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp) / spec_summ
								end do

								do dim1 = 1, dimensions
									constructor%v_f(dim1,dim,i,j,k) =	constructor%v%v_ptr%pr(dim1)%cells(i,j,k)
								end do

							end if
					!	end if
				end do										
				end do	
				end do	
			end do

			do spec = 1,species_number
				constructor%Y_f_new%v_ptr%pr(spec)%cells(:,:,:,:) = constructor%Y_f(spec,:,:,:,:)			
			end do
			
			do dim = 1,dimensions
				constructor%v_f_new%v_ptr%pr(dim)%cells(:,:,:,:) = constructor%v_f(dim,:,:,:,:)    
			end do
		
			constructor%p_f_new%s_ptr%cells(:,:,:,:)	= constructor%p_f     
			constructor%rho_f_new%s_ptr%cells(:,:,:,:)	= constructor%rho_f   
			constructor%E_f_f_new%s_ptr%cells(:,:,:,:)	= constructor%E_f_f	
			constructor%v_s_f_new%s_ptr%cells(:,:,:,:)	= constructor%v_s_f
			
		end if

		call constructor%mpi_support%exchange_flow_scalar_field(constructor%p_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%rho_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%E_f_f_new%s_ptr)	
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%v_s_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_vector_field(constructor%Y_f_new%v_ptr)
		call constructor%mpi_support%exchange_flow_vector_field(constructor%v_f_new%v_ptr)
	
		do spec = 1,species_number
			constructor%Y_f(spec,:,:,:,:)	= constructor%Y_f_new%v_ptr%pr(spec)%cells		
		end do
		
		do dim = 1,dimensions
			constructor%v_f(dim,:,:,:,:)     = constructor%v_f_new%v_ptr%pr(dim)%cells
		end do
	

		constructor%time		        =   calculation_time
		constructor%time_step	        =   manager%solver_options%get_initial_time_step()
		constructor%initial_time_step   =   manager%solver_options%get_initial_time_step()
		
		call constructor%state_eq%apply_state_equation_flow_variables_for_IC()

		call constructor%apply_boundary_conditions_main()
		
        constructor%p_f     = constructor%p_f_new%s_ptr%cells
        constructor%rho_f   = constructor%rho_f_new%s_ptr%cells
		constructor%E_f_f	= constructor%E_f_f_new%s_ptr%cells
		constructor%v_s_f	= constructor%v_s_f_new%s_ptr%cells

        allocate(flame_front_coords(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))
        
        open(newunit = flame_loc_unit, file = 'av_flame_data.dat', status = 'replace', form = 'formatted')
        
call manager%create_timer(cabaret_timer                 ,'CABARET solver time'              , 'sol_t')
        call manager%create_timer(cabaret_gas_dynamics_timer    ,'CABARET gas dynamics time'        , 'gd_t')
        call manager%create_timer(cabaret_eos_timer             ,'CABARET eos solver time'          , 'eos_t')
        call manager%create_timer(cabaret_chemistry_timer       ,'CABARET chemistry solver time'    , 'chem_t')
        call manager%create_timer(cabaret_diffusion_timer       ,'CABARET diffusion solver time'    , 'diff_t')
        call manager%create_timer(cabaret_heattransfer_timer    ,'CABARET heattransfer solver time' , 'ht_t')
        call manager%create_timer(cabaret_viscosity_timer       ,'CABARET viscosity solver time'    , 'visc_t')
	end function

	subroutine solve_test_problem(this)
		class(cabaret_solver)	,intent(inout)	:: this

		associate(	rho			=> this%rho%s_ptr, &
					rho_f_new	=> this%rho_f_new%s_ptr, &
					v_f_new		=> this%v_f_new%v_ptr, &
					p			=> this%p%s_ptr	) 

			rho%cells = this%domain%get_processor_rank()
			rho_f_new%cells = this%domain%get_processor_rank()
			v_f_new%pr(1)%cells	= this%domain%get_processor_rank()
			v_f_new%pr(2)%cells	= this%domain%get_processor_rank()*10

			call this%mpi_support%exchange_conservative_scalar_field(rho)

			call this%mpi_support%exchange_flow_scalar_field(rho_f_new)

			call this%mpi_support%exchange_flow_vector_field(v_f_new)

		end associate

	end subroutine


	subroutine solve_problem(this)
		class(cabaret_solver)	,intent(inout)	:: this
    
		real(dp)	,dimension(2)							:: r_inv, q_inv, r_inv_corrected, q_inv_corrected, r_inv_new, q_inv_new
		real(dp)	,dimension(:,:)	,allocatable	,save	:: v_inv, v_inv_corrected, v_inv_new
		real(dp)	,dimension(:)	,allocatable	,save	:: v_inv_half, v_inv_old
		real(dp)	,dimension(:,:)	,allocatable	,save	:: Y_inv, Y_inv_corrected, y_inv_new
		real(dp)	,dimension(:)	,allocatable	,save	:: Y_inv_half, Y_inv_old

		real(dp)					:: r_inv_half, q_inv_half, R_inv_old, Q_inv_old
		real(dp)					:: G_half, G_half_old	, G_half_lower, G_half_higher
		
		real(dp)	,dimension(3)	:: characteristic_speed

		real(dp)	:: v_f_approx, v_s_f_approx
		real(dp)	:: v_f_approx_lower, v_f_approx_higher

        real(dp)	:: r
		real(dp)	:: g_inv, alpha = 0.0_dp, alpha_loc
		real(dp)	:: f, corr, diss_l, diss_r, diss = 0.5_dp
		integer		,save :: dissipator_active = 0
		real(dp) :: max_inv, min_inv, maxmin_inv
		real(dp) :: mean_higher, mean_lower
		real(dp)	:: sources
		real(dp)	:: mean_sources
		real(dp)	:: summ_frac
		real(dp)	:: energy_output_time	= 2.0e-07_dp
		real(dp)	:: energy_output_radii	= 4.0e-04_dp
		real(dp)	:: energy_source = 1.9e+04_dp
		real(dp)	,save	:: energy_output_rho = 0.0_dp
		real(dp)	,save	:: energy_output = 0.0_dp
		integer		,save	:: energy_output_flag = 0 
		
		!real(dp)	:: r_inf, q_inf, s_inf
		real(dp)	,parameter	:: u_inf		= 0.0
		real(dp)	,parameter	:: p_inf		= 101325.000000000_dp
		real(dp)	,parameter	:: rho_inf		= 0.487471044003493_dp
		real(dp)	,parameter	:: c_inf		= 539.709011600784_dp
		real(dp)	,parameter	:: gamma_inf	= 1.40137_dp
		real(dp)			:: g_inf		= 1.0_dp / rho_inf / c_inf
		real(dp)	,save	:: q_inf		= u_inf - p_inf / rho_inf / c_inf
	
		real(dp)	:: spec_summ, rho_Y
		
		real(dp), save	:: Mach_old
		real(dp)			:: Mach
		real(dp)			:: flow_Mach, cons_Mach
		real(dp)			:: Mach_rise, Max_Mach_rise
		integer				:: Max_Mach_rise_i, Max_Mach_rise_j
		real(dp)			:: Max_v_s, Min_v_s
		
		integer				:: nu
		
		integer             :: sign, bound_number

		integer	:: dimensions, species_number
		character(len=20)	:: coordinate_system
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		real(dp)	,dimension(3)	:: cell_size
        integer     	:: particles_phase_counter

		integer			:: i,j,k,plus,dim,dim1,dim2,spec,iter		

		integer	:: thread
        
        logical	:: stabilized
		
		call cabaret_timer%tic()
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number

		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()

		flow_utter_loop = this%domain%get_local_utter_faces_bounds()
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		if (.not.allocated(v_inv)) then
			allocate(v_inv(dimensions,2), v_inv_corrected(dimensions,2), v_inv_new(dimensions,2))
			allocate(v_inv_half(dimensions), v_inv_old(dimensions))
			allocate(Y_inv(species_number,2), Y_inv_corrected(species_number,2), y_inv_new(species_number,2))
			allocate(Y_inv_half(species_number), Y_inv_old(species_number))
		end if
	
		call this%apply_boundary_conditions_main()						

		select case(coordinate_system)
        case ('cartesian')	
				nu = 1
			case ('cylindrical')
				! x -> z, y -> r
				nu = 2
			case ('spherical')
				! x -> r
				nu = 3 
		end select		
		
        associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &

					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old)

		call this%mpi_support%exchange_conservative_scalar_field(p)
		call this%mpi_support%exchange_conservative_scalar_field(rho)
		call this%mpi_support%exchange_conservative_scalar_field(v_s)
		call this%mpi_support%exchange_conservative_scalar_field(E_f)

		call this%mpi_support%exchange_conservative_vector_field(Y)
		call this%mpi_support%exchange_conservative_vector_field(v)

		call cabaret_gas_dynamics_timer%tic()

		rho_old		= rho%cells
		E_f_old		= E_f%cells
		v_s_old		= v_s%cells
		p_old		= p%cells

		do spec = 1,species_number
			Y_old(spec,:,:,:)		=	Y%pr(spec)%cells
		end do

		do dim = 1,dimensions
			v_old(dim,:,:,:)		=	v%pr(dim)%cells
		end do	

		end associate


		associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					gamma		=> this%gamma%s_ptr		, &

                    v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &

                    v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					p_f			=> this%p_f		, &
					Y_f			=> this%Y_f		, &
					
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	,&

					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

        !$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r)
		!!$omp& shared(this,cons_inner_loop,dimensions,species_number,Max_v_s,Min_v_s,cell_size,nu,coordinate_system)

		!$omp do collapse(3) schedule(guided) reduction(max:Max_v_s) reduction(min:Min_v_s)	 		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				
                r = mesh%mesh(1,i,j,k)
                
				rho%cells(i,j,k)	= 0.0_dp
				E_f%cells(i,j,k)	= 0.0_dp
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k)	=	0.0_dp 
				end do
		  
				do dim = 1,dimensions
					v%pr(dim)%cells(i,j,k)	=	0.0_dp 		
                end do				

				do dim = 1,dimensions
                    
                    mean_higher	= rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) 
					mean_lower	= rho_f(dim,i,j,k) * v_f(dim,dim,i,j,k)                   

					rho%cells(i,j,k)	=	rho%cells(i,j,k) - (mean_higher - mean_lower) /cell_size(1)
                    
                    if(dim == 1) then
						rho%cells(i,j,k)	=	rho%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																	*	0.5_dp * (mean_higher +	mean_lower)																		&
																	*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																	/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 	
                    end if
				end do
				
				rho%cells(i,j,k)		=	rho_old(i,j,k)	+  0.5_dp * this%time_step * rho%cells(i,j,k)

				spec_summ = 0.0_dp
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k)		= 0.0_dp
					do	dim = 1,dimensions
						mean_higher	=  rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * y_f(spec,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f(dim,i,j,k) * y_f(spec,dim,i,j,k) * v_f(dim,dim,i,j,k)
						
						Y%pr(spec)%cells(i,j,k)	=  Y%pr(spec)%cells(i,j,k)	-	(mean_higher - mean_lower ) /cell_size(1)
                        
                        if(dim == 1) then
							Y%pr(spec)%cells(i,j,k)	=  Y%pr(spec)%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																				*	0.5_dp * (mean_higher +	mean_lower)																		&
																				*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																				/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 	
                        end if
					end do
                    
					Y%pr(spec)%cells(i,j,k)		=	rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * Y%pr(spec)%cells(i,j,k)
					Y%pr(spec)%cells(i,j,k)		=	Y%pr(spec)%cells(i,j,k)	/ rho%cells(i,j,k)
					
					spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dp)
				end do                
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dp) / spec_summ 
				end do

				do dim = 1,dimensions
					do dim1 = 1,dimensions
                        
						mean_higher	=  rho_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))*v_f(dim,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))*v_f(dim1,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))
						mean_lower	=  rho_f(dim1,i,j,k)*v_f(dim,dim1,i,j,k)*v_f(dim1,dim1,i,j,k)
                        
						v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower ) /cell_size(1)

                        if(dim1 == 1) then                        
							v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																				*	0.5_dp * (mean_higher +	mean_lower)																		&
																				*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																				/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 
                        end if                        

					end do

					v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k)	-	 ( p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p_f(dim,i,j,k)) /cell_size(1)
					
					v%pr(dim)%cells(i,j,k)	=	rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dp*this%time_step*v%pr(dim)%cells(i,j,k) 
					v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k) /rho%cells(i,j,k)
	
				end do	
    
				do dim = 1,dimensions
                    
 					mean_higher	=  (rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*E_f_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+	p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))*v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
					mean_lower	=  (rho_f(dim,i,j,k)*E_f_f(dim,i,j,k)																	+	p_f(dim,i,j,k))									*v_f(dim,dim,i,j,k)                   
                    
					E_f%cells(i,j,k)	= 	E_f%cells(i,j,k) -	(mean_higher - mean_lower ) /cell_size(1)
                    
                    if(dim == 1) then
						E_f%cells(i,j,k)	=	E_f%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																	*	0.5_dp * (mean_higher +	mean_lower)																		&
																	*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																	/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 
                    end if
				end do	
				
				E_f%cells(i,j,k)		=	rho_old(i,j,k) * E_f_old(i,j,k)  + 0.5_dp*this%time_step* E_f%cells(i,j,k)
				E_f%cells(i,j,k)		=	E_f%cells(i,j,k) / rho%cells(i,j,k) 
				
			end if
		end do
		end do
        end do
		!$omp end do nowait		
		!$omp end parallel				
        
                    end associate
        
!        call this%check_symmetry()
		
        call cabaret_gas_dynamics_timer%toc()
		! ********** Conservative variables state eq *******************	
        call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation() 
        call cabaret_eos_timer%toc(new_iter=.true.)
		! **************************************************************
        call cabaret_gas_dynamics_timer%tic()
        
        associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					gamma		=> this%gamma%s_ptr		, &
     
					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr)

		! Block one mpi_exchange
		call this%mpi_support%exchange_conservative_scalar_field(p)
		call this%mpi_support%exchange_conservative_scalar_field(rho)
		call this%mpi_support%exchange_conservative_scalar_field(E_f)
		call this%mpi_support%exchange_conservative_scalar_field(v_s)
		call this%mpi_support%exchange_conservative_scalar_field(gamma)

		
		call this%mpi_support%exchange_conservative_vector_field(Y)
		call this%mpi_support%exchange_conservative_vector_field(v)

		end associate

		associate(	p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr)

				!$omp parallel default(shared)  private(i,j,k,dim,dim1,spec) !, &
		!!$omp& shared(this,flow_utter_loop,species_number,dimensions)

		do dim = 1,dimensions
			!$omp do collapse(3) schedule(guided)			
			do k = flow_utter_loop(3,1),flow_utter_loop(3,2)
			do j = flow_utter_loop(2,1),flow_utter_loop(2,2)
			do i = flow_utter_loop(1,1),flow_utter_loop(1,2)		
			
				p_f_new%cells(dim,i,j,k)	= 0.0_dp
				rho_f_new%cells(dim,i,j,k)	= 0.0_dp	
		
				do dim1 = 1,dimensions
					v_f_new%pr(dim1)%cells(dim,i,j,k)	= 0.0_dp
				end do
				do spec = 1,species_number
					Y_f_new%pr(spec)%cells(dim,i,j,k)	= 0.0_dp		
				end do		
			end do	
			end do
			end do
			!$omp end do nowait
		end do
		!$omp end parallel		
		
        end associate		
        
        associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					v_s			=> this%v_s%s_ptr		, &

					v			=> this%v%v_ptr	, &

                    v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					p_f			=> this%p_f		, &
					Y_f			=> this%Y_f		, &
					
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	,&
					
					r_inv_corr      => this%r_inv_corr	, &
					q_inv_corr      => this%q_inv_corr	, &
					v_inv_corr      => this%v_inv_corr	,&            

					bc				=> this%boundary%bc_ptr)
					
		!$omp parallel default(shared)  private(thread,i,j,k,dim,dim1,loop,G_half,G_half_old,G_half_lower,G_half_higher,r_inv,R_inv_half,R_inv_old,q_inv,Q_inv_half,Q_inv_old,v_inv,v_inv_half,v_inv_old,r_inv_new,q_inv_new,v_inv_new,g_inv,max_inv,min_inv,maxmin_inv,v_f_approx,v_s_f_approx,characteristic_speed,diss_l,diss_r,alpha_loc,sign,bound_number) !, &
		!!$omp& shared(this,cons_utter_loop,cons_inner_loop,dimensions,Max_v_s,Min_v_s,diss,alpha,dissipator_active,cell_size,lock,coordinate_system)
		
		do dim = 1,dimensions

			thread = 0
		
			loop(3,1) = cons_utter_loop(3,1)*I_m(dim,3) + cons_inner_loop(3,1)*(1 - I_m(dim,3))
			loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

			loop(2,1) = cons_utter_loop(2,1)*I_m(dim,2) + cons_inner_loop(2,1)*(1 - I_m(dim,2))
			loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

			loop(1,1) = cons_utter_loop(1,1)*I_m(dim,1) + cons_inner_loop(1,1)*(1 - I_m(dim,1))
			loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))						
			
			!$omp do collapse(3) schedule(guided)				
			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)

				if(bc%bc_markers(i,j,k) == 0) then

					G_half			= 1.0_dp / (v_s%cells(i,j,k)*rho%cells(i,j,k))
					G_half_old		= 1.0_dp / (v_s_old(i,j,k)*rho_old(i,j,k))

					! *********** Riemann quasi invariants *************************

					! ********* Lower invariants ***********
					
					r_inv(1) 	= v_f(dim,dim,i,j,k)									+ G_half*p_f(dim,i,j,k)
					r_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
					R_inv_half	= v%pr(dim)%cells(i,j,k)								+ G_half*(p%cells(i,j,k))
					R_inv_old	= v_old(dim,i,j,k)										+ G_half_old*(p_old(i,j,k))
					
					q_inv(1) 	= v_f(dim,dim,i,j,k)									- G_half*p_f(dim,i,j,k)
					q_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
					Q_inv_half	= v%pr(dim)%cells(i,j,k)								- G_half*(p%cells(i,j,k))
					Q_inv_old	= v_old(dim,i,j,k)										- G_half_old*(p_old(i,j,k))

					do dim1 = 1,dimensions
						if (dim1 == dim) then
							v_inv(dim,1)		= p_f(dim,i,j,k)                                    - v_s%cells(i,j,k)**2*rho_f(dim,i,j,k)
							v_inv(dim,2)		= p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- v_s%cells(i,j,k)**2*rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 
							v_inv_half(dim)		= p%cells(i,j,k)                                    - v_s%cells(i,j,k)**2*rho%cells(i,j,k)
							v_inv_old(dim)		= p_old(i,j,k)                                      - v_s_old(i,j,k)**2*rho_old(i,j,k)
						else
							v_inv(dim1,1)		= v_f(dim1,dim,i,j,k)
							v_inv(dim1,2)		= v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
							v_inv_half(dim1)	= v%pr(dim1)%cells(i,j,k)
							v_inv_old(dim1)		= v_old(dim1,i,j,k)
						end if
					end do

					! ******************* Linear interpolation *********************

					diss_l = 0.00_dp
					diss_r = 0.00_dp

					r_inv_new(1) = (2.0_dp*R_inv_half - (1.0_dp-diss_l)*r_inv(2))/(1.0_dp+diss_l)
					q_inv_new(1) = (2.0_dp*Q_inv_half - (1.0_dp-diss_l)*q_inv(2))/(1.0_dp+diss_l)
	
					r_inv_new(2) = (2.0_dp*R_inv_half - (1.0_dp-diss_r)*r_inv(1))/(1.0_dp+diss_r)
					q_inv_new(2) = (2.0_dp*Q_inv_half - (1.0_dp-diss_r)*q_inv(1))/(1.0_dp+diss_r)				

					do dim1 = 1,dimensions
						v_inv_new(dim1,1)	= (2.0_dp*v_inv_half(dim1) - (1.0_dp-diss_l)*v_inv(dim1,2))/(1.0_dp+diss_l)
						v_inv_new(dim1,2)	= (2.0_dp*v_inv_half(dim1) - (1.0_dp-diss_r)*v_inv(dim1,1))/(1.0_dp+diss_r)
					end do

					! **************** Non-linear flow correction ******************

					alpha_loc = 0.0_dp

					g_inv = 0.0_dp
					g_inv = ((R_inv_half - R_inv_old)/(0.5_dp*this%time_step) + (v%pr(dim)%cells(i,j,k) + v_s%cells(i,j,k))*(r_inv(2) - r_inv(1))/cell_size(1))
				!	g_inv = ((R_half - R_old)/(this%time_step) + (v%pr(dim)%cells(i,j,k) + v_s%cells(i,j,k))*(r(2) - r(1))/cell_size(1))

					max_inv = max(r_inv(1),R_inv_half,r_inv(2)) + g_inv*this%time_step
					min_inv = min(r_inv(1),R_inv_half,r_inv(2)) + g_inv*this%time_step
					maxmin_inv = abs(max_inv - min_inv)

					max_inv = max_inv + (-alpha_loc)*maxmin_inv
					min_inv = min_inv - (-alpha_loc)*maxmin_inv

					if ((min_inv <= r_inv_new(1)).and.(r_inv_new(1) <= max_inv))    r_inv_corr(1,dim,i,j,k) = r_inv_new(1)
					if (r_inv_new(1) < min_inv)										r_inv_corr(1,dim,i,j,k) = min_inv
					if (max_inv < r_inv_new(1))										r_inv_corr(1,dim,i,j,k) = max_inv
					
					if ((min_inv <= r_inv_new(2)).and.(r_inv_new(2) <= max_inv))    r_inv_corr(2,dim,i,j,k) = r_inv_new(2)
					if (r_inv_new(2) < min_inv)										r_inv_corr(2,dim,i,j,k) = min_inv
					if (max_inv < r_inv_new(2))										r_inv_corr(2,dim,i,j,k) = max_inv               
					
					g_inv = 0.0_dp
					g_inv = ((Q_inv_half - Q_inv_old)/(0.5_dp*this%time_step) + (v%pr(dim)%cells(i,j,k) - v_s%cells(i,j,k))*(q_inv(2) - q_inv(1))/cell_size(1))
                !   g_inv = ((Q_half - Q_old)/(this%time_step) + (v%pr(dim)%cells(i,j,k) - v_s%cells(i,j,k))*(q(2) - q(1))/cell_size(1))

					max_inv = max(q_inv(1),Q_inv_half,q_inv(2)) + g_inv*this%time_step
					min_inv = min(q_inv(1),Q_inv_half,q_inv(2)) + g_inv*this%time_step
					maxmin_inv = abs(max_inv - min_inv)
				
					max_inv = max_inv + (-alpha_loc)*maxmin_inv
					min_inv = min_inv - (-alpha_loc)*maxmin_inv
	
					if ((min_inv <= q_inv_new(1)).and.(q_inv_new(1) <= max_inv))    q_inv_corr(1,dim,i,j,k) = q_inv_new(1)
					if (q_inv_new(1) < min_inv)										q_inv_corr(1,dim,i,j,k) = min_inv
					if (max_inv < q_inv_new(1))										q_inv_corr(1,dim,i,j,k) = max_inv
						
					if ((min_inv <= q_inv_new(2)).and.(q_inv_new(2) <= max_inv))	q_inv_corr(2,dim,i,j,k) = q_inv_new(2)
					if (q_inv_new(2) < min_inv)										q_inv_corr(2,dim,i,j,k) = min_inv
					if (max_inv < q_inv_new(2))										q_inv_corr(2,dim,i,j,k) = max_inv

					do dim2 = 1,dimensions
						g_inv = 0.0_dp

						g_inv = ((v_inv_half(dim2) - v_inv_old(dim2))/(0.5_dp*this%time_step) + (v%pr(dim)%cells(i,j,k))*(v_inv(dim2,2) - v_inv(dim2,1))/cell_size(1))
                    !   g_inv = ((v_inv_half(dim2) - v_inv_old(dim2))/(this%time_step) + (v%pr(dim)%cells(i,j,k))*(v_inv(dim2,2) - v_inv(dim2,1))/cell_size(1))
                        if(dim2 == dim) g_inv = 0.0_dp
                        
                        
						max_inv = max(v_inv(dim2,1),v_inv_half(dim2),v_inv(dim2,2)) + g_inv*this%time_step
						min_inv = min(v_inv(dim2,1),v_inv_half(dim2),v_inv(dim2,2)) + g_inv*this%time_step
						maxmin_inv = abs(max_inv - min_inv)

						max_inv = max_inv + (-alpha_loc)*maxmin_inv
						min_inv = min_inv - (-alpha_loc)*maxmin_inv
						
						if ((min_inv <= v_inv_new(dim2,1)).and.(v_inv_new(dim2,1) <= max_inv))		v_inv_corr(1,dim2,dim,i,j,k) = v_inv_new(dim2,1)
						if (v_inv_new(dim2,1) < min_inv)											v_inv_corr(1,dim2,dim,i,j,k) = min_inv
						if (max_inv < v_inv_new(dim2,1))											v_inv_corr(1,dim2,dim,i,j,k) = max_inv
					
						if ((min_inv <= v_inv_new(dim2,2)).and.(v_inv_new(dim2,2) <= max_inv))		v_inv_corr(2,dim2,dim,i,j,k) = v_inv_new(dim2,2)
						if (v_inv_new(dim2,2) < min_inv)											v_inv_corr(2,dim2,dim,i,j,k) = min_inv
						if (max_inv < v_inv_new(dim2,2))											v_inv_corr(2,dim2,dim,i,j,k) = max_inv   											
		
                    end do
                    
                    !**************************** Boundary conditions *****************************
					if (( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,1) ) .and. &
						( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,2) )) then
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then

								v_f_approx		= 0.5_dp * (v%pr(dim)%cells(i,j,k)	+ v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
								v_s_f_approx	= 0.5_dp * (v_s%cells(i,j,k)			+ v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))								
		
								characteristic_speed(1) = v_f_approx + v_s_f_approx
								characteristic_speed(2) = v_f_approx - v_s_f_approx
								characteristic_speed(3) = v_f_approx		
									
								call this%apply_boundary_conditions_flow(dim, i,j,k, characteristic_speed, q_inv_corr(:,dim,i,j,k), r_inv_corr(:,dim,i,j,k), v_inv_corr(:,:,dim,i,j,k), G_half)
								
							end if
						end do
					end if
				end if
			end do
			end do
			end do
			!$omp end do
		
        end do

		!$omp end parallel                    
        end associate
                    
        associate(	rho			=> this%rho%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
            
					v			=> this%v%v_ptr	, &
					
					v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					p_f			=> this%p_f		, &
					
					r_inv_corr      => this%r_inv_corr	, &
					q_inv_corr      => this%q_inv_corr	, &
					v_inv_corr      => this%v_inv_corr	,&            

					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &
					
					bc			=> this%boundary%bc_ptr)
            
		!$omp parallel default(shared)  private(thread,i,j,k,dim,dim1,loop,G_half_lower,G_half_higher,v_f_approx,v_s_f_approx,characteristic_speed,sign,bound_number) !, &
		!!$omp& shared(this,cons_utter_loop,flow_inner_loop,dimensions,cell_size,lock,coordinate_system)        

        do dim = 1, dimensions

			loop = flow_inner_loop

			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))	
			end do		

		!$omp do collapse(3) schedule(guided)

			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)
				if ((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then        

					G_half_lower	= 1.0_dp / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))*rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					G_half_higher	= 1.0_dp / (v_s%cells(i,j,k)*rho%cells(i,j,k))
                        
					v_f_approx		= 0.5_dp*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v%pr(dim)%cells(i,j,k))
					v_s_f_approx	= 0.5_dp*(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))			+ v_s%cells(i,j,k))                    
                    
					if	((.not.((abs(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	<	abs(	v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))	.and.(	  v%pr(dim)%cells(i,j,k)	>		v_s%cells(i,j,k))))	&
					.and.(.not.((	 v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	<			-v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	.and.(abs(v%pr(dim)%cells(i,j,k))	< abs(	v_s%cells(i,j,k)))))) then 

						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx
					
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) >= 0.0_dp )) then				

							p_f_new%cells(dim,i,j,k)			=	(r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	-   q_inv_corr(1,dim,i,j,k)) / (G_half_lower + G_half_higher)
							rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k) - v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) / v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(G_half_higher * r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + G_half_lower * q_inv_corr(1,dim,i,j,k))	/ (G_half_lower + G_half_higher)
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
  								end if
                            end do
                            
                            if(characteristic_speed(3) == 0.0_dp) then
                                rho_f_new%cells(dim,i,j,k)		= rho_f(dim,i,j,k)
                            end if
						end if		
			
						if (( characteristic_speed(1) >= 0.0_dp ).and.&
							( characteristic_speed(2) < 0.0_dp ).and.&
							( characteristic_speed(3) < 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)			=	(r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	-   q_inv_corr(1,dim,i,j,k)) / (G_half_lower + G_half_higher)
							rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k) - v_inv_corr(1,dim,dim,i,j,k)) / (v_s%cells(i,j,k)**2)
							
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(G_half_higher * r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + G_half_lower * q_inv_corr(1,dim,i,j,k))	/ (G_half_lower + G_half_higher)
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
								end if
							end do
							continue
						end if

						if (( characteristic_speed(1) >= 0.0_dp ).and.&
							( characteristic_speed(2) >= 0.0_dp ).and.&
							( characteristic_speed(3) >= 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)	= 0.5_dp * (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - q_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / G_half_lower
							rho_f_new%cells(dim,i,j,k)	= (p_f_new%cells(dim,i,j,k) - v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) / v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
						
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	0.5_dp * (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + q_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
									!if ( characteristic_speed(3) == 0 ) then
									!	v_f_new%pr(dim1)%cells(dim,i,j,k) = v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									!end if
								end if
							end do	
						end if
					
						if (( characteristic_speed(1) < 0.0_dp ).and.&
							( characteristic_speed(2) < 0.0_dp ).and.&
							( characteristic_speed(3) < 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)			= 0.5_dp * (r_inv_corr(1,dim,i,j,k) - q_inv_corr(1,dim,i,j,k)) / G_half_higher
							rho_f_new%cells(dim,i,j,k)			= (p_f_new%cells(dim,i,j,k) - v_inv_corr(1,dim,dim,i,j,k)) / v_s%cells(i,j,k)**2
						
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	0.5_dp * (r_inv_corr(1,dim,i,j,k) + q_inv_corr(1,dim,i,j,k))
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
								end if
							end do
                        end if
                            
						if (( characteristic_speed(1) == 0.0_dp ).and.&
							( characteristic_speed(2) == 0.0_dp ).and.&
							( characteristic_speed(3) == 0.0_dp )) then  
                            

                            p_f_new%cells(dim,i,j,k)		= p_f(dim,i,j,k)
                            rho_f_new%cells(dim,i,j,k)		= rho_f(dim,i,j,k)
                            
                            do dim1 = 1,dimensions
								v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_f(dim1,dim,i,j,k)
							end do
                        end if
                            
                    end if
                    
					!**************************** Sound points *****************************
     
					v_f_approx		= 0.5_dp*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v%pr(dim)%cells(i,j,k))
					v_s_f_approx	= 0.5_dp*(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))			+ v_s%cells(i,j,k))
     
					characteristic_speed(1) = v_f_approx + v_s_f_approx
					characteristic_speed(2) = v_f_approx - v_s_f_approx
					characteristic_speed(3) = v_f_approx
     
					if (((abs(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	<	abs(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))) &
					.and.(	  v%pr(dim)%cells(i,j,k)	>	v_s%cells(i,j,k))) )then 
						
                            
						v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.5_dp*(v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) &
																*0.5_dp*(v_s%cells(i,j,k) + v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
		   
						p_f_new%cells(dim,i,j,k)			= (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - v_f_new%pr(dim)%cells(dim,i,j,k))/G_half_lower
     
                            
						if (characteristic_speed(3) >= 0.0_dp) then
							rho_f_new%cells(dim,i,j,k)		= (p_f_new%cells(dim,i,j,k) - v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2)
							
							do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= (p_f_new%cells(dim,i,j,k) - v_inv_corr(1,dim,dim,i,j,k)) / (v_s%cells(i,j,k)**2)
                                
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                                end if
							end do
                        end if                            
						
                        continue
                            
					end if
     
					if (((	 v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	 <	  -v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	&
					.and.	(abs(v%pr(dim)%cells(i,j,k)) < abs(v_s%cells(i,j,k)))))then 
					
						v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.5_dp*(v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) &
																*0.5_dp*(v_s%cells(i,j,k) + v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
		   
						p_f_new%cells(dim,i,j,k)			= (v_f_new%pr(dim)%cells(dim,i,j,k) - q_inv_corr(1,dim,i,j,k))/G_half_higher
						
						if (characteristic_speed(3) >= 0.0_dp) then
							rho_f_new%cells(dim,i,j,k)		= (p_f_new%cells(dim,i,j,k) - v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2)
							
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= (p_f_new%cells(dim,i,j,k) - v_inv_corr(1,dim,dim,i,j,k)) / (v_s%cells(i,j,k)**2)
                                
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                                end if
							end do
                                
                        end if   
                            
                        continue
					end if
                end if
			end do
			end do
			end do
			!$omp end do
		
		end do
		!$omp end parallel		
		
        end associate

        ! ************************************************  		

                    
        associate(  p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr)

		call this%mpi_support%exchange_flow_scalar_field(rho_f_new)
		call this%mpi_support%exchange_flow_scalar_field(p_f_new)

        end associate

        associate(	v			=> this%v%v_ptr	            , &
					Y			=> this%Y%v_ptr	            , &
                    Y_f			=> this%Y_f		            , &
					Y_old		=> this%Y_old	            , &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					bc			=> this%boundary%bc_ptr)
        
		!$omp parallel default(shared)  private(i,j,k,dim,loop,spec,Y_inv,Y_inv_half,Y_inv_new,Y_inv_old,g_inv,max_inv,min_inv,maxmin_inv,alpha_loc,Y_inv_corrected,v_f_approx_lower,v_f_approx_higher,spec_summ,bound_number,diss_r,diss_l) !, &
		!!$omp& shared(this,cons_utter_loop,cons_inner_loop,dimensions,species_number,Max_v_s,Min_v_s,diss,alpha,cell_size,lock)
		
		do dim = 1,dimensions

			! Avoid looping in transverse direction in ghost cells

			loop(3,1) = cons_utter_loop(3,1)*I_m(dim,3) + cons_inner_loop(3,1)*(1 - I_m(dim,3))
			loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

			loop(2,1) = cons_utter_loop(2,1)*I_m(dim,2) + cons_inner_loop(2,1)*(1 - I_m(dim,2))
			loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

			loop(1,1) = cons_utter_loop(1,1)*I_m(dim,1) + cons_inner_loop(1,1)*(1 - I_m(dim,1))
			loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))						

			!$omp do collapse(3) schedule(guided)	
			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)		
			
				if(bc%bc_markers(i,j,k) == 0) then
	
					do spec = 1,species_number
						y_inv(spec,1)		= y_f(spec,dim,i,j,k)
						y_inv(spec,2)		= y_f(spec,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						y_inv_half(spec)	= y%pr(spec)%cells(i,j,k)
						y_inv_old(spec)		= y_old(spec,i,j,k)
                        
                      
						!y_inv(spec,1)		= Y_f(spec,dim,i,j,k)									* rho%cells(i,j,k)		- Y%pr(spec)%cells(i,j,k) / v_s%cells(i,j,k)**2 * p_f(dim,i,j,k)							
						!y_inv(spec,2)		= Y_f(spec,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* rho%cells(i,j,k)		- Y%pr(spec)%cells(i,j,k) / v_s%cells(i,j,k)**2 * p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
						!Y_inv_half(spec)	= Y%pr(spec)%cells(i,j,k)								* rho%cells(i,j,k)		- Y%pr(spec)%cells(i,j,k) / v_s%cells(i,j,k)**2 * p%cells(i,j,k)		
						!Y_inv_old(spec)		= Y_old(spec,i,j,k)										* rho_old(i,j,k)		- Y%pr(spec)%cells(i,j,k) / v_s%cells(i,j,k)**2 * p_old(i,j,k)	                        

                    end do
					
					diss_l = 0.0_dp
					diss_r = 0.0_dp

					do spec = 1,species_number
						y_inv_new(spec,1)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_l)*y_inv(spec,2))/(1.0_dp+diss_l)
						y_inv_new(spec,2)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_r)*y_inv(spec,1))/(1.0_dp+diss_r)
					end do	
					
					do spec = 1,species_number
						g_inv =  0.0_dp
						
						g_inv = ((Y_inv_half(spec) - Y_inv_old(spec))/(0.5_dp*this%time_step) + (v%pr(dim)%cells(i,j,k))*(y_inv(spec,2) - y_inv(spec,1))/cell_size(1))
					!	g_inv = ((Y_inv_half(spec) - Y_inv_old(spec))/(this%time_step) + (v%pr(dim)%cells(i,j,k))*(y_inv(spec,2) - y_inv(spec,1))/cell_size(1))
						
                        
						alpha_loc = 0.0_dp
						
						max_inv	= max(y_inv(spec,1),Y_inv_half(spec),y_inv(spec,2)) +  g_inv*this%time_step
						min_inv	= min(y_inv(spec,1),Y_inv_half(spec),y_inv(spec,2)) +  g_inv*this%time_step

						maxmin_inv = abs(max_inv - min_inv)
						
						max_inv = max_inv + (-alpha_loc)*maxmin_inv
						min_inv = min_inv - (-alpha_loc)*maxmin_inv

                        
						max_inv	= min(max_inv,1.0_dp)
						min_inv	= max(min_inv,0.0_dp)

						if ((min_inv <= y_inv_new(spec,1)).and.(y_inv_new(spec,1) <= max_inv))		y_inv_corrected(spec,1) = y_inv_new(spec,1)
						if (y_inv_new(spec,1) < min_inv)											y_inv_corrected(spec,1) = min_inv
						if (max_inv < y_inv_new(spec,1))											y_inv_corrected(spec,1) = max_inv
					
						if ((min_inv <= y_inv_new(spec,2)).and.(y_inv_new(spec,2) <= max_inv))		y_inv_corrected(spec,2) = y_inv_new(spec,2)
						if (y_inv_new(spec,2) < min_inv)											y_inv_corrected(spec,2) = min_inv
						if (max_inv < y_inv_new(spec,2))											y_inv_corrected(spec,2) = max_inv
						
					end do					
					
#ifdef OMP					
					call omp_set_lock(lock(i,j,k))
					call omp_set_lock(lock(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))	
#endif

					!# Boundary conditions on Y are set up along with other flow variables, here values on boundaries should not be updated
					do spec = 1,species_number
						if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,1) ) then
							bound_number	= bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
							if(bound_number == 0) then	
							!	v_f_approx_lower		= v_f_new%pr(dim)%cells(dim,i,j,k) 
                                v_f_approx_lower		= 0.5_dp * (v%pr(dim)%cells(i,j,k)	+ v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								if (v_f_approx_lower < 0.0_dp) then
                                !	Y_f_new%pr(spec)%cells(dim,i,j,k) =  (y_inv_corrected(spec,1) + Y%pr(spec)%cells(i,j,k) * p_f(dim,i,j,k)  / v_s%cells(i,j,k)**2 )  / rho%cells(i,j,k)
                                	Y_f_new%pr(spec)%cells(dim,i,j,k) =  (y_inv_corrected(spec,1))
								end if	
							end if
						end if

						if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,2) ) then
							bound_number	= bc%bc_markers(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
							if(bound_number == 0) then
								!v_f_approx_higher		= v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                                v_f_approx_higher		= 0.5_dp * (v%pr(dim)%cells(i,j,k)	+ v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
								if (v_f_approx_higher >= 0.0_dp) then
                                !	Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (y_inv_corrected(spec,2) + Y%pr(spec)%cells(i,j,k) * p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  / v_s%cells(i,j,k)**2 )  / rho%cells(i,j,k)
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (y_inv_corrected(spec,2)) 	
   								end if
                            end if
                        end if
					end do
					
#ifdef OMP						
					call omp_unset_lock(lock(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
					call omp_unset_lock(lock(i,j,k))	
#endif
				end if
			end do
			end do
			end do
		!$omp end do
			
		!$omp do collapse(3) schedule(guided)	
			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)	
				if(bc%bc_markers(i,j,k) == 0) then
					if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,1) ) 	then
						spec_summ = 0.0_dp
						do spec = 1,species_number
							spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
						end do

						do spec = 1,species_number
						!	Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp) / spec_summ
						end do
					end if
				end if
			end do
			end do
			end do
			!$omp end do nowait	

        end do
		!$omp end parallel

        end associate
                    
!        call this%check_symmetry()        

                    
        associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					
		            v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
                
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	)
            
        rho_old		= rho%cells
		E_f_old		= E_f%cells
		v_s_old		= v_s%cells
		p_old		= p%cells

		do spec = 1,species_number
			Y_old(spec,:,:,:)		=	Y%pr(spec)%cells
		end do

		do dim = 1,dimensions
			v_old(dim,:,:,:)		=	v%pr(dim)%cells
		end do	        
        
                    end associate
        
        call cabaret_gas_dynamics_timer%toc()
		! ******************* Eqn of state ***************
        call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation_flow_variables() 
        call cabaret_eos_timer%toc(new_iter=.true.)		
        ! ************************************************  
        call cabaret_gas_dynamics_timer%tic()

        
        associate(	v_s_f_new	=> this%v_s_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					e_i_f_new	=> this%e_i_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_scalar_field(e_i_f_new)
		call this%mpi_support%exchange_flow_scalar_field(E_f_f_new)	
		call this%mpi_support%exchange_flow_scalar_field(v_s_f_new)
		call this%mpi_support%exchange_flow_vector_field(v_f_new)		
		call this%mpi_support%exchange_flow_vector_field(Y_f_new)

                    end associate

        ! *********** Conservative variables calculation ***************
		
		associate(	rho			=> this%rho%s_ptr		, &
					E_f			=> this%E_f%s_ptr		, &

		        	v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&

					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &

					bc				=> this%boundary%bc_ptr     , &
					mesh			=> this%mesh%mesh_ptr)
            
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r) !, &
		!!$omp& shared(this,cons_inner_loop,dimensions,species_number,cell_size,nu,coordinate_system)
         
		!$omp do collapse(3) schedule(guided)			
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
  			
			if (bc%bc_markers(i,j,k) == 0) then	
                
                r = mesh%mesh(1,i,j,k)
                
				rho%cells(i,j,k)	= 0.0_dp
				E_f%cells(i,j,k)	= 0.0_dp
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k)	=	0.0_dp 
				end do
		  
				do dim = 1,dimensions
					v%pr(dim)%cells(i,j,k)	=	0.0_dp 		
                end do	
                
				do dim = 1,dimensions
                    
                    mean_higher	= rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) 
					mean_lower	= rho_f_new%cells(dim,i,j,k)    *v_f_new%pr(dim)%cells(dim,i,j,k)                    

					rho%cells(i,j,k)	=	rho%cells(i,j,k) - (mean_higher - mean_lower) /cell_size(1)
                    
                    if(dim == 1) then
						rho%cells(i,j,k)	=	rho%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																	*	0.5_dp * (mean_higher +	mean_lower)																		&
																	*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																	/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1)))
                    end if
				end do                

				rho%cells(i,j,k)	=	rho_old(i,j,k) + 0.5_dp*this%time_step*rho%cells(i,j,k)
	   
				spec_summ = 0.0_dp
				do spec = 1,species_number
					do	dim = 1,dimensions
						mean_higher	=  rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  * y_f_new%pr(spec)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f_new%cells(dim,i,j,k)    * y_f_new%pr(spec)%cells(dim,i,j,k)   *v_f_new%pr(dim)%cells(dim,i,j,k)
						
						Y%pr(spec)%cells(i,j,k)	=  Y%pr(spec)%cells(i,j,k)	-	(mean_higher - mean_lower ) /cell_size(1)
                        
                        if(dim == 1) then
							Y%pr(spec)%cells(i,j,k)	=  Y%pr(spec)%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																				*	0.5_dp * (mean_higher +	mean_lower)																		&
																				*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																				/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 
                        end if
					end do
					
					Y%pr(spec)%cells(i,j,k)		=	rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * Y%pr(spec)%cells(i,j,k)
					Y%pr(spec)%cells(i,j,k)		=	Y%pr(spec)%cells(i,j,k)	/ rho%cells(i,j,k)
					
					spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dp)
				end do
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dp) / spec_summ 
				end do

				do dim = 1,dimensions
					do dim1 = 1,dimensions
						mean_higher	= rho_f_new%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))  *v_f_new%pr(dim)%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) *v_f_new%pr(dim1)%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))
						mean_lower	= rho_f_new%cells(dim1,i,j,k)  *v_f_new%pr(dim)%cells(dim1,i,j,k) *v_f_new%pr(dim1)%cells(dim1,i,j,k)
													
						v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower)	/cell_size(1)
                        
                        
                        if(dim1 == 1) then                        
							v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																				*	0.5_dp * (mean_higher +	mean_lower)																		&
																				*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																				/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 
                        end if                        
					end do
	   
					mean_higher	= p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= p_f_new%cells(dim,i,j,k)
													
					v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower)	/cell_size(1)
					
					v%pr(dim)%cells(i,j,k)	= rho_old(i,j,k)*v_old(dim,i,j,k) +  0.5_dp * this%time_step * v%pr(dim)%cells(i,j,k)
					v%pr(dim)%cells(i,j,k)	= v%pr(dim)%cells(i,j,k) / rho%cells(i,j,k) 
				end do	
	   
				E_f%cells(i,j,k)		=	0.0_dp  
				do dim = 1,dimensions
					mean_higher	= (rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))*E_f_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))	+	p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))*v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= (rho_f_new%cells(dim,i,j,k)*E_f_f_new%cells(dim,i,j,k)	+	p_f_new%cells(dim,i,j,k))*v_f_new%pr(dim)%cells(dim,i,j,k)	
				
					E_f%cells(i,j,k)        = 	E_f%cells(i,j,k)		-	(mean_higher - mean_lower)	/cell_size(1)
                    
                    if(dim == 1) then
						E_f%cells(i,j,k)        = 	E_f%cells(i,j,k)		-	2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																				*	0.5_dp * (mean_higher +	mean_lower)																		&
																				*	((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))							&
																				/	((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1))) 
                    end if
				end do	
				E_f%cells(i,j,k) = rho_old(i,j,k) * E_f_old(i,j,k) + 0.5_dp * this%time_step * E_f%cells(i,j,k)
				E_f%cells(i,j,k) = E_f%cells(i,j,k) /rho%cells(i,j,k)
	
			end if
        end do
		end do
		end do
        ! **************************************************************
		!$omp end do nowait		
		!$omp end parallel

        end associate

!        call this%check_symmetry()

        associate(	v   => this%v%v_ptr)
			call this%mpi_support%exchange_conservative_vector_field(v)	
        end associate

		call this%apply_boundary_conditions_main()
        
        call cabaret_gas_dynamics_timer%toc()

		if(this%additional_particles_phases_number /= 0) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
    			call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_main(this%time)
			end do	            
        end if
                
        call cabaret_heattransfer_timer%tic()
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)
        call cabaret_heattransfer_timer%toc(new_iter=.true.)
        
        call cabaret_diffusion_timer%tic()
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
        call cabaret_diffusion_timer%toc(new_iter=.true.)
        
        call cabaret_viscosity_timer%tic()
		if (this%viscosity_flag)	call this%viscosity_solver%solve_viscosity(this%time_step)
        call cabaret_viscosity_timer%toc(new_iter=.true.)
        
        call cabaret_chemistry_timer%tic()
		if (this%reactive_flag)		call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
        call cabaret_chemistry_timer%toc(new_iter=.true.)
        

		if(this%additional_particles_phases_number /= 0) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
				call this%particles_solver(particles_phase_counter)%particles_solve(this%time_step)				!# Lagrangian particles solver
!				call this%particles_solver(particles_phase_counter)%particles_euler_step_v_E(this%time_step)	!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%particles_lagrange_step(this%time_step)		!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%particles_final_step(this%time_step)		!# Continuum particles solver		
			end do		
		end if 	        
        
        call cabaret_gas_dynamics_timer%tic()
        
        associate(  E_f_prod		=> this%E_f_prod			, &
			        v_prod			=> this%v_prod				, &
			        Y_prod			=> this%Y_prod				, &
			        rho_prod		=> this%rho_prod)
    
		E_f_prod	= 0.0_dp
		rho_prod	= 0.0_dp
		Y_prod		= 0.0_dp
		v_prod		= 0.0_dp	
        end associate
        
        associate(	rho			=> this%rho%s_ptr		, &
					E_f			=> this%E_f%s_ptr		, &

		        	v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &

                    v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff	=> this%E_f_prod_diff%s_ptr	, &
					E_f_prod		=> this%E_f_prod			, &
					v_prod			=> this%v_prod				, &
					Y_prod			=> this%Y_prod				, &
					rho_prod		=> this%rho_prod			, &
            
       
					rho_prod_particles  => this%rho_prod_particles	, &
					E_f_prod_particles  => this%E_f_prod_particles	, &
                    v_prod_particles	=> this%v_prod_particles	, &
					Y_prod_particles	=> this%Y_prod_particles	, &

					bc				=> this%boundary%bc_ptr         , &
					mesh			=> this%mesh%mesh_ptr)
		
		!$omp parallel default(shared)  private(i,j,k,dim,spec,spec_summ) !, &
		!!$omp& shared(this,cons_inner_loop,species_number,dimensions,energy_source,energy_output_rho,energy_output,energy_output_flag,energy_output_time,energy_output_radii,cell_size,coordinate_system)
    		
		!$omp do collapse(3) schedule(guided)		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			
			if(bc%bc_markers(i,j,k) == 0) then	
				
				if (this%reactive_flag)	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_chem%cells(i,j,k) * this%time_step
					do spec = 1,species_number
						Y_prod(spec,i,j,k)	= Y_prod(spec,i,j,k)	+ Y_prod_chem%pr(spec)%cells(i,j,k) * this%time_step
					end do		
                end if
                
				if (this%heat_trans_flag)	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_heat%cells(i,j,k) * this%time_step
                end if
                
 				if (this%diffusion_flag)	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_diff%cells(i,j,k) * this%time_step
					do spec = 1, species_number
						Y_prod(spec,i,j,k)	= Y_prod(spec,i,j,k)	+ Y_prod_diff%pr(spec)%cells(i,j,k) * this%time_step
					end do
				end if

				if (this%viscosity_flag)	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_visc%cells(i,j,k) * this%time_step
					do dim = 1, dimensions
						v_prod(dim,i,j,k)	= v_prod(dim,i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k) * this%time_step
				!		v_prod(dim,i,j,k)	= v_prod(dim,i,j,k) + this%g(dim) * (rho%cells(1,1,1) - rho%cells(i,j,k)) * this%time_step
					end do
                end if		
                
				! ************************* Energy release ******************
				if (energy_output_flag == 1) then 
					if(this%time <= energy_output_time) then
						if(mesh%mesh(1,i,j,k) <= mesh%mesh(1,1,j,k) + energy_output_radii) then
							E_f_prod(i,j,k)			= E_f_prod(i,j,k)	+ energy_source * 1.0e+10 * this%time_step *  rho%cells(i,j,k) !* cell_size(1) ! * 4.0_dp * Pi * mesh%mesh(1,i,j,k) * mesh%mesh(1,i,j,k)
							energy_output_rho		= energy_output_rho	+ energy_source * 1.0e+10 * this%time_step *  rho%cells(i,j,k) * cell_size(1) * 4.0_dp * Pi  ! * mesh%mesh(1,i,j,k) * mesh%mesh(1,i,j,k)
							energy_output			= energy_output		+ energy_source * 1.0e+10 * this%time_step
						end if
					else	
						energy_output_flag = 2
					end if
				else
					if (energy_output_flag == 2) then
						print *, ' Energy input	: ', energy_output_rho
						print *, ' Time	: ', this%time
						print *, 'r_0 : ', mesh%mesh(1,1,j,k), ' r_f : ', energy_output_radii/cell_size(1)
						stop
						energy_output_flag = 0
					end if				
				end if
				! ***********************************************************				
				
				E_f%cells(i,j,k) = E_f%cells(i,j,k) + E_f_prod(i,j,k)/rho%cells(i,j,k)
				
				spec_summ = 0.0_dp
				do spec = 1, species_number
					Y%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k) + Y_prod(spec,i,j,k)/rho%cells(i,j,k)                   
					spec_summ = spec_summ + Y%pr(spec)%cells(i,j,k)
                end do		
				
				do dim = 1, dimensions
					v%pr(dim)%cells(i,j,k)	= v%pr(dim)%cells(i,j,k) + v_prod(dim,i,j,k) /rho%cells(i,j,k)
                end do					
                
                if (this%additional_particles_phases_number /= 0) then
                    spec_summ = 0.0_dp
					do particles_phase_counter = 1, this%additional_particles_phases_number
						E_f%cells(i,j,k)	= E_f%cells(i,j,k) + E_f_prod_particles(particles_phase_counter)%s_ptr%cells(i,j,k)
                        rho%cells(i,j,k)    = rho%cells(i,j,k) + rho_prod_particles(particles_phase_counter)%s_ptr%cells(i,j,k)
                        do dim = 1, dimensions
                            v%pr(dim)%cells(i,j,k)	= v%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k)
                        end do
                        do spec = 1, species_number
							Y%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k) + Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k)
                            spec_summ = spec_summ + Y%pr(spec)%cells(i,j,k)
						end do	
					end do		
                end if                 
                                
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dp) / spec_summ 
				end do
                
			end if	
		end do	
		end do
		end do
		!$omp end do nowait
		!$omp end parallel	
				
        end associate
        
        call cabaret_gas_dynamics_timer%toc()

        call cabaret_eos_timer%tic()		
		call this%state_eq%apply_state_equation() 		
		call cabaret_eos_timer%toc(new_iter=.true.)
        
        call cabaret_gas_dynamics_timer%tic()
		associate(  v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					p_f			=> this%p_f		, &
					v_s_f		=> this%v_s_f	, & 
                    Y_f			=> this%Y_f		, &
     !       
					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					v_s_f_new	=> this%v_s_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr)
            
		!$omp parallel default(shared)  private(i,j,k,dim,spec) !, &
		!!$omp& shared(this,flow_utter_loop,species_number,dimensions)

		do dim = 1,dimensions
			!$omp do collapse(3) schedule(static)			
			do k = flow_utter_loop(3,1),flow_utter_loop(3,2)
			do j = flow_utter_loop(2,1),flow_utter_loop(2,2)
			do i = flow_utter_loop(1,1),flow_utter_loop(1,2)		
		
				do spec = 1,species_number
					Y_f(spec,dim,i,j,k)	= Y_f_new%pr(spec)%cells(dim,i,j,k)		
				end do

				do dim1 = 1,dimensions
					v_f(dim1,dim,i,j,k)		= v_f_new%pr(dim1)%cells(dim,i,j,k)
				end do
	
				p_f(dim,i,j,k)	    = p_f_new%cells(dim,i,j,k)	
				rho_f(dim,i,j,k)	= rho_f_new%cells(dim,i,j,k)	
				E_f_f(dim,i,j,k)	= E_f_f_new%cells(dim,i,j,k)
				v_s_f(dim,i,j,k)	= v_s_f_new%cells(dim,i,j,k)

			end do
			end do
			end do
			! **************************************************************
			!$omp end do nowait	
        end do
		!$omp end parallel
		
        end associate

		!call this%state_eq%check_conservation_laws()

        !call this%if_stabilized(this%time, stabilized)
        !if (stabilized) stop		

        if (this%CFL_condition_flag) then
			call this%calculate_time_step()
        end if

		this%time	= this%time + this%time_step
        
        call cabaret_gas_dynamics_timer%toc(new_iter=.true.)
        call cabaret_timer%toc(new_iter=.true.)

	end subroutine

	subroutine apply_boundary_conditions_flow(this, dim,i,j,k, characteristic_speed, q_inv_corrected, r_inv_corrected, v_inv_corrected, G_half)

		class(cabaret_solver)		,intent(inout)		:: this
		integer						,intent(in)			:: i, j, k, dim
		real(dp)	,dimension(3)	,intent(in)			:: characteristic_speed
		real(dp)	,dimension(2)	,intent(in)			:: q_inv_corrected, r_inv_corrected
		real(dp)	,dimension(:,:)	,intent(in)			:: v_inv_corrected
		real(dp)					,intent(in)			:: G_half
		
		real(dp)				:: r_inf, q_inf, s_inf, G_half_inf

		real(dp)				:: spec_summ

		real(dp)	,dimension(3)	:: cell_size
		character(len=20)		:: boundary_type_name
		integer					:: dimensions, species_number
		character(len=20)		:: coordinate_system
		integer					:: sign, bound_number
		integer 				:: plus, dim1, spec	

		associate(  v_s				=> this%v_s%s_ptr			, &
					rho				=> this%rho%s_ptr			, &
					p				=> this%p%s_ptr				, &
					E_f				=> this%E_f%s_ptr			, &
					v_f				=> this%v_f					, &
					v				=> this%v%v_ptr				, &
					Y				=> this%Y%v_ptr				, &
					p_f_new			=> this%p_f_new%s_ptr		, &	
					rho_f_new		=> this%rho_f_new%s_ptr		, &
					E_f_f_new		=> this%E_f_f_new%s_ptr		, &
					v_f_new			=> this%v_f_new%v_ptr		, &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		
		do plus = 1,2
			sign			= (-1)**plus
			bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
			if( bound_number /= 0 ) then
				boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
				select case(boundary_type_name)
					case('wall','symmetry_plane')
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) >= 0.0_dp )) then	
							!# Subsonic flow from left to right, wall on the left.
							if (sign == -1) then
								v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.0_dp 
								do dim1 = 1,dimensions
									if( dim1 /= dim) then							
										v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v%pr(dim1)%cells(i,j,k)
									end if
								end do		
								p_f_new%cells(dim,i,j,k)			=	-q_inv_corrected(1)	/ G_half
								rho_f_new%cells(dim,i,j,k)			=	rho%cells(i,j,k)  !/ mesh%mesh(1,i,j,k) / mesh%mesh(1,i,j,k)  * (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)) * (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))
								
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k)	= 	Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp) / spec_summ
								end do
							end if
							!# Subsonic flow from left to right, wall on the right.
							if (sign == 1) then
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_inv_corrected(2,dim1)
									end if
								end do
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= r_inv_corrected(2)/G_half
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(2,dim)) / (v_s%cells(i,j,k)**2)
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do								
							end if
						end if
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) < 0.0_dp )) then
							!# Subsonic flow from right to left, wall on the left.
							if (sign == -1) then
								p_f_new%cells(dim,i,j,k)			=	-q_inv_corrected(1)/ G_half
								rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k)	-	v_inv_corrected(1,dim))	/ (v_s%cells(i,j,k)**2)
								v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.0_dp
								do dim1 = 1,dimensions
									if( dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corrected(1,dim1)
									end if
								end do
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k)	= 	Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp) / spec_summ
								end do
							end if
							!# Subsonic flow from right to left, wall on the right.
							if (sign == 1) then
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v%pr(dim1)%cells(i,j,k)
									end if
								end do
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= rho%cells(i,j,k)
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= r_inv_corrected(2)/G_half
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do	

							end if

						end if
									
					case ('outlet')
								
						! ******************* Acoustic outlet ***********************************
					
						G_half_inf = 1.0_dp / (rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
						r_inf = v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + G_half_inf * p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	
						q_inf = v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - G_half_inf * p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	
						s_inf = p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2 * rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
					
						
						if (( characteristic_speed(1) >= 0.0_dp )		.and.&
							( characteristic_speed(2) < 0.0_dp ))then	!.and.&
							if (sign == -1) then											!# Lower edge (left,bottom,far), flow from lower to higher edge  
								
								p_f_new%cells(dim,i,j,k)			=	(r_inf - q_inv_corrected(1))/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i,j,k)	=	(G_half*r_inf + G_half_inf*q_inv_corrected(1))/(G_half_inf + G_half)
								rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k) - s_inf)/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2		
								
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i,j,k)	= 0.0_dp
									end if
								end do
								
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k)	= 	Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp) / spec_summ
								end do
								
								E_f_f_new%cells(dim,i,j,k)			=	E_f%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								
							end if
							if (sign == 1) then
							
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(r_inv_corrected(2) - q_inf)/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	(G_half_inf*r_inv_corrected(2) + G_half*q_inf)/(G_half_inf + G_half)
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(2,dim))/v_s%cells(i,j,k)**2								
							
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_inv_corrected(2,dim1)
									end if
								end do

								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do								
							end if
						end if
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) < 0.0_dp )) then

							if (sign == 1) then
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(r_inv_corrected(2) - q_inf)/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	(G_half_inf*r_inv_corrected(2) + G_half*q_inf)/(G_half_inf + G_half)
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - s_inf)/v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))**2								
							
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.0_dp
									end if
								end do

								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do	

							end if

						end if						

												
						! ******************* Shock outlet ***********************************
							
						if (( characteristic_speed(1) >= 0.0_dp ).and.&
							( characteristic_speed(2) >= 0.0_dp ).and.&
							( characteristic_speed(3) >= 0.0_dp )) then
											
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.5_dp * (r_inv_corrected(2) - q_inv_corrected(2)) / G_half
											
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(2,dim)) / (v_s%cells(i,j,k)**2)
											
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.5_dp * (r_inv_corrected(2) + q_inv_corrected(2))
								else
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_inv_corrected(2,dim1)
									if ( characteristic_speed(3) == 0 ) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									end if
								end if
							end do
						end if

					case ('inlet')
							
						G_half_inf = 1.0_dp / (rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))
						r_inf = v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) + G_half_inf * p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	
						q_inf = v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - G_half_inf * p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	
						s_inf = p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2 * rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
					
						! ******************* Acoustic inlet ***********************************					
						if (( characteristic_speed(1) >= 0.0_dp )		.and.&
							( characteristic_speed(2) < 0.0_dp ))then	!.and.&
							if (sign == -1) then											!# Lower edge (left,bottom,far), flow from lower to higher edge  
								
								p_f_new%cells(dim,i,j,k)			=	(r_inf - q_inv_corrected(1))/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i,j,k)	=	(G_half*r_inf + G_half_inf*q_inv_corrected(1))/(G_half_inf + G_half)
								rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k) - s_inf)/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2		
								
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i,j,k)	= 0.0_dp
									end if
								end do
								
								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k)	= 	Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp) / spec_summ
								end do
								
								E_f_f_new%cells(dim,i,j,k)			=	E_f%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								
							end if
							if (sign == 1) then
							
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(r_inv_corrected(2) - q_inf)/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	(G_half_inf*r_inv_corrected(2) + G_half*q_inf)/(G_half_inf + G_half)
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(2,dim))/v_s%cells(i,j,k)**2								
							
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_inv_corrected(2,dim1)
									end if
								end do

								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do								
							end if
						end if
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) < 0.0_dp )) then
							if (sign == 1) then
								p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(r_inv_corrected(2) - q_inf)/(G_half_inf + G_half)
								v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	(G_half_inf*r_inv_corrected(2) + G_half*q_inf)/(G_half_inf + G_half)
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - s_inf)/v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))**2								
							
								do dim1 = 1, dimensions
									if (dim1 /= dim) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.0_dp
									end if
								end do

								spec_summ = 0.0_dp
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 	= Y%pr(spec)%cells(i,j,k)  	
									spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
								end do
								do spec = 1,species_number
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
								end do	

							end if

						end if					
																		
						! ******************* Shock inlet ***********************************
							
						if (( characteristic_speed(1) >= 0.0_dp ).and.&
							( characteristic_speed(2) >= 0.0_dp ).and.&
							( characteristic_speed(3) >= 0.0_dp )) then
											
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.5_dp * (r_inv_corrected(2) - q_inv_corrected(2)) / G_half
											
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(2,dim)) / (v_s%cells(i,j,k)**2)
											
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.5_dp * (r_inv_corrected(2) + q_inv_corrected(2))
								else
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_inv_corrected(2,dim1)
									if ( characteristic_speed(3) == 0 ) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									end if
								end if
                            end do
                            
                            spec_summ = 0.0_dp
							do spec = 1,species_number
								Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=  Y%pr(spec)%cells(i,j,k) 	
								spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp)
							end do
							do spec = 1,species_number
								Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = max(Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), 0.0_dp) / spec_summ
							end do
						end if
	
				end select
			end if
		end do
		end associate

	end subroutine

	subroutine apply_boundary_conditions_main(this)

		class(cabaret_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_utter_loop, cons_inner_loop
		character(len=20)		:: boundary_type_name
		real(dp)				:: farfield_density, farfield_pressure, wall_temperature, farfield_velocity
        real(dp)             :: relax_time, delay_time, H_wall_frac

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,dim2,specie_number

		dimensions			= this%domain%get_domain_dimensions()
		cons_utter_loop		= this%domain%get_local_utter_cells_bounds()	
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()	
		
        H_wall_frac			= 0.0_dp
        
		associate(  T				=> this%T%s_ptr					, &
					mol_mix_conc	=> this%mol_mix_conc%s_ptr		, &
					p				=> this%p%s_ptr					, &
					rho				=> this%rho%s_ptr				, &
					v				=> this%v%v_ptr					, &
					v_f_new			=> this%v_f_new%v_ptr			, &
					v_s				=> this%v_s%s_ptr				, &
					Y				=> this%Y%v_ptr					, &
					bc				=> this%boundary%bc_ptr			, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(shared)  private(i,j,k,plus,dim,dim1,sign,bound_number,farfield_pressure,farfield_density,farfield_velocity, relax_time, delay_time, wall_temperature,boundary_type_name) !, &
		!!$omp& shared(this,dimensions,cons_utter_loop,cons_inner_loop)
            
		!$omp do collapse(3) schedule(static)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							!if(((i+sign)*I_m(dim,1) + (j+sign)*I_m(dim,2) + (k+sign)*I_m(dim,3) <= cons_utter_loop(dim,2)).and. &
							!   ((i+sign)*I_m(dim,1) + (j+sign)*I_m(dim,2) + (k+sign)*I_m(dim,3) >= cons_utter_loop(dim,1))) then

								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then

									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')

											p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= p%cells(i,j,k)
											rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i,j,k)
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= T%cells(i,j,k)
											mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= mol_mix_conc%cells(i,j,k)

											v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= v_s%cells(i,j,k)
									
											do dim1 = 1, dimensions
												if(dim1 == dim) then
													v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v%pr(dim1)%cells(i,j,k)
												else
													v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k)
												end if
											end do

											do specie_number = 1, this%chem%chem_ptr%species_number
												Y%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(specie_number)%cells(i,j,k)
											end do
			
											if(bc%boundary_types(bound_number)%is_conductive()) then
												wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
											end if
											if(.not.bc%boundary_types(bound_number)%is_slip()) then
												do dim1 = 1, dimensions
													if(dim1 /= dim) then
														v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - 1.0_dp * v%pr(dim1)%cells(i,j,k) 
														! v%pr(dim1)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - 6.0_dp * v%pr(dim1)%cells(i,j,k) 
														!- 10.0_dp * v%pr(dim1)%cells(i,j,k)
													end if
												end do
											end if


										case ('outlet')
										!	farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										!	farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()
										!	v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  sign*sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
                                        case ('inlet')
                                            !if(bound_number == 5) then
                                            !    delay_time = 50e-06_dp
                                            !    relax_time = 20e-06_dp
                                            !    farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
                                            !    v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	farfield_velocity*(min((this%time-delay_time)/relax_time,1.0))
                                            !end if
										!	farfield_pressure	= bc%boundary_types(bound_number)%get_farfield_pressure()
										!	farfield_density	= bc%boundary_types(bound_number)%get_farfield_density()
										!	v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -sign*sqrt(abs((p%cells(i,j,k) - farfield_pressure)*(rho%cells(i,j,k) - farfield_density)/farfield_density/rho%cells(i,j,k)))
									end select

								end if
							!end if
						end do
					end do
				end if
			end do
			end do
			end do

		!$omp end do nowait
		!$omp end parallel

		end associate
	end subroutine
	
	subroutine calculate_time_step(this)

#ifdef mpi
	use MPI
#endif

		class(cabaret_solver)	,intent(inout)	:: this
		
		real(dp)	:: delta_t_interm, time_step(1), velocity_value
		real(dp)	,dimension(:)	,allocatable	,save	:: time_step_array

		integer						:: dimensions
		integer						:: processor_rank, processor_number, mpi_communicator

		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		integer	:: sign
		integer :: i,j,k,dim,error

		processor_rank		= this%domain%get_processor_rank()
		mpi_communicator	= this%domain%get_mpi_communicator()

		if (.not.allocated(time_step_array)) then
			processor_number = this%domain%get_mpi_communicator_size()
			allocate(time_step_array(processor_number))
            print *, 'time_step_array allocation'
		end if

		time_step(1)	= 1.0e-05_dp !this%initial_time_step

		associate(  v				=> this%v%v_ptr		, &
					v_s				=> this%v_s%s_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)
		
		dimensions			= this%domain%get_domain_dimensions()
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()
		cell_size			= mesh%get_cell_edges_length()					
					
		!!$omp parallel default(shared)  private(i,j,k,dim,delta_t_interm,velocity_value) , &
		!!$omp& firstprivate(this)	,&
		!!$omp& shared(v,v_s,mesh,bc,time_step)
		!!$omp do collapse(3) schedule(static) reduction(min:time_step)
					
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dp
				do dim = 1,dimensions
					velocity_value = velocity_value + v%pr(dim)%cells(i,j,k)*v%pr(dim)%cells(i,j,k)
				end do
				delta_t_interm = minval(cell_size,cell_size > 0.0_dp) / (sqrt(velocity_value) + v_s%cells(i,j,k))
				if (delta_t_interm < time_step(1)) then
					time_step(1) = delta_t_interm
				end if
			end if
		end do
		end do
		end do
	
		!!$omp end do nowait
		!!$omp end parallel

		time_step_array(processor_rank+1) = time_step(1) 
		
#ifdef mpi					
		call mpi_gather(time_step,1,MPI_DOUBLE_PRECISION,time_step_array,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif
		
		if (processor_rank == 0) then
			do i = 0,size(time_step_array) - 1 
				if (time_step_array(i+1) < time_step(1)) time_step(1) = time_step_array(i+1)
			end do
		end if

#ifdef mpi	
		call mpi_bcast(time_step,1,MPI_DOUBLE_PRECISION,0,mpi_communicator,error)
#endif

		this%time_step = this%courant_fraction * time_step(1)

		if(time_step(1) < 5.7e-09_dp) then
			this%time_step = 0.25 * time_step(1)
			print *, 'Time step was reduced. Co = 1: ',time_step(1), '. Reduced: ', this%time_step
		end if

		end associate
			
	end subroutine
	
	subroutine set_CFL_coefficient(this,coefficient)
		class(cabaret_solver)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: coefficient
	
		this%courant_fraction = coefficient
		
	end subroutine
	
	pure function get_time_step(this)
		real(dp)						:: get_time_step
		class(cabaret_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dp)						:: get_time
		class(cabaret_solver)	,intent(in)		:: this

		get_time = this%time
	end function

	subroutine check_symmetry(this)
    
		class(cabaret_solver)	,intent(inout)	:: this
		
		real(dp)	:: delta_t_interm, time_step(1), velocity_value
		real(dp)	,dimension(:)	,allocatable	,save	:: time_step_array

		integer						:: dimensions, species_number
		integer						:: processor_rank, processor_number, mpi_communicator

		integer		,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
        integer		,dimension(3,2)	:: loop
		real(dp)	,dimension(3)	:: cell_size
		integer	:: sign
		integer :: i,j,k,dim,dim1,error, spec

		associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v			=> this%v%v_ptr			, &
					Y			=> this%Y%v_ptr			, &
					
					v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					e_i_f		=> this%e_i_f	, &
					p_f			=> this%p_f		, &
					v_s_f		=> this%v_s_f	, & 
					Y_f			=> this%Y_f		, &
					
					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &

					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)
		
		dimensions			= this%domain%get_domain_dimensions()
		species_number      = this%chem%chem_ptr%species_number
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()
        cons_utter_loop		= this%domain%get_local_utter_cells_bounds()
        flow_inner_loop		= this%domain%get_local_inner_faces_bounds()

		cell_size			= mesh%get_cell_edges_length()					
					
        !# central symmetry
        
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				if (rho%cells(i,j,k) - rho%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
                    print *, "density asymmetry", i,j,k, rho%cells(i,j,k) - rho%cells(cons_inner_loop(1,2)-i,cons_inner_loop(2,2)-j,cons_inner_loop(3,2)-k)
                    pause
                end if
                if (E_f%cells(i,j,k) - E_f%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
                    print *, "energy asymmetry", i,j,k, E_f%cells(i,j,k) - E_f%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1)
                    pause
                end if
                do dim = 1, dimensions
					if (abs(v%pr(dim)%cells(i,j,k)) - abs(v%pr(dim)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1)) /= 0 ) then
						print *, "velocity ",dim, " asymmetry", i,j,k, abs(v%pr(dim)%cells(i,j,k)) - abs(v%pr(dim)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1))
                        pause
                    end if 
				end do
                do spec = 1, species_number
					if (Y%pr(spec)%cells(i,j,k) - Y%pr(spec)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
						print *, "Concentration ",spec, " asymmetry", i,j,k, Y%pr(spec)%cells(i,j,k) - Y%pr(spec)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) 
                        pause
                    end if 
				end do            
			end if
		end do
		end do
        end do

		do dim = 1,dimensions

			! Avoid looping in transverse direction in ghost cells

			loop(3,1) = cons_inner_loop(3,1)*I_m(dim,3) + cons_inner_loop(3,1)*(1 - I_m(dim,3))
			loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

			loop(2,1) = cons_inner_loop(2,1)*I_m(dim,2) + cons_inner_loop(2,1)*(1 - I_m(dim,2))
			loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

			loop(1,1) = cons_inner_loop(1,1)*I_m(dim,1) + cons_inner_loop(1,1)*(1 - I_m(dim,1))
			loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))						

			!$omp do collapse(3) schedule(guided)	
			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2)
			do i = loop(1,1),loop(1,2)
                    
				if (rho_f_new%cells(dim,i,j,k) - rho_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
					print *, "flow density asymmetry" 
                    print *, i,j,k, loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1
                    print *, rho_f_new%cells(dim,i,j,k) - rho_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
					pause
				end if
				if (p_f_new%cells(dim,i,j,k) - p_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
					print *, "flow pressure asymmetry", i,j,k, p_f_new%cells(dim,i,j,k) - p_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
					pause
				end if
				do dim1 = 1, dimensions
					if (abs(v_f_new%pr(dim)%cells(dim,i,j,k)) - abs(v_f_new%pr(dim)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)) /= 0 ) then
						print *, "Flow velocity ",dim, " asymmetry", i,j,k, abs(v_f_new%pr(dim)%cells(dim,i,j,k)) - abs(v_f_new%pr(dim)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1))
						pause
					end if 
				end do
				do spec = 1, species_number
					if (Y_f_new%pr(spec)%cells(dim,i,j,k) - Y_f_new%pr(spec)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
						print *, "Flow concentration ",spec, " asymmetry", i,j,k, Y_f_new%pr(spec)%cells(dim,i,j,k) - Y_f_new%pr(spec)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
						pause
					end if 
                end do

			end do
			end do
			end do        
        
        end do
        
		end associate    
    
    
    
    end subroutine
    
    subroutine if_stabilized(this,time,stabilized)
		class(cabaret_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)			:: time    

		logical				,intent(out)	:: stabilized
		
        logical						:: boundary 
		real(dp)	,dimension(3)	:: cell_size		
		
		real(dp)					:: max_val, left_val, right_val, flame_velocity, flame_surface_length, surface_factor
		real(dp)					:: a, b 
		real(dp)					:: time_diff, time_delay, time_stabilization
		real(dp), save			:: previous_flame_location = 0.0_dp, current_flame_location = 0.0_dp, farfield_velocity = 0.0_dp
		real(dp), save			:: previous_time = 0.0_dp, current_time = 0.0_dp
        real(dp), save			:: av_flame_velocity = 0.0_dp, previous_av_flame_velocity = 0.0_dp
		real(dp), dimension(20),	save	:: flame_velocity_array = 0.0_dp
		integer		,save			:: correction = 0, counter = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dp)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2 
		
		integer :: CO_index, H2O2_index, HO2_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		
		stabilized = .false.
		
		associate (	v				=> this%v%v_ptr				, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr)
	
		time_delay			= 1e-05_dp			
		time_diff			= 1e-04_dp
		time_stabilization	= 5e-06_dp			
		
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# 1D front tracer
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				max_val = 0.0_dp
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
					
						!! Grad temp
						!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val) then
						!	max_val = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
						!	flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
						!	flame_front_index = i
						!end if
					
						! max HO2
						if (abs(Y%pr(HO2_index)%cells(i,j,k)) > max_val) then
							max_val = Y%pr(HO2_index)%cells(i,j,k)
							flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
							flame_front_index = i
						end if	
					end if
				end do
			end do
			end do
	
			!left_val	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_val	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
			 
			left_val	= Y%pr(HO2_index)%cells(flame_front_index-1,1,1)
			right_val	= Y%pr(HO2_index)%cells(flame_front_index+1,1,1)				 
			
			a = (right_val + left_val - 2.0_dp * max_val)/2.0_dp/cell_size(1)**2
			b = (max_val - left_val)/cell_size(1) - a*(2.0_dp*flame_front_coords(1) - cell_size(1))
			
			current_flame_location = -b/2.0_dp/a

			if(correction == 0) then
				previous_flame_location = current_flame_location
            end if
			
            boundary = .false.
            if(flame_front_index > cons_inner_loop(1,2) - 10) boundary = .true.
            
			if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
                
				flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)

				previous_flame_location = current_flame_location
				previous_time = current_time
				
                previous_av_flame_velocity = sum(flame_velocity_array)
                
                do i = 19, 1, -1
					flame_velocity_array(i+1) = flame_velocity_array(i)  
                end do
                flame_velocity_array(1) = flame_velocity
                
                av_flame_velocity = sum(flame_velocity_array)
                
				if( (correction /= 0).and.(abs(av_flame_velocity - previous_av_flame_velocity) < 1e-03))then 
					counter = counter + 1
                end if
                
                write (flame_loc_unit,'(5E14.6)') time, current_flame_location, flame_velocity, av_flame_velocity, abs(av_flame_velocity - previous_av_flame_velocity)

			end if

			if ((counter > 10).or.boundary) then
				stabilized = .true.
			end if
			
			correction = correction + 1
			
		end if	
			
		end associate
    end subroutine	
    
    
end module
