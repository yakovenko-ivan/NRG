module fds_low_mach_solver_class

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
	use chemical_kinetics_solver_class
	use chemical_properties_class
	
	use viscosity_solver_class
	use fourier_heat_transfer_solver_class
	use fickean_diffusion_solver_class
	use lagrangian_droplets_solver_class
	use lagrangian_particles_solver_class
	use droplets_solver_class
	
	use solver_options_class
	
	implicit none
	
#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: fds_solver, fds_solver_c

	type(field_scalar_cons)	,target	:: p_dyn	,p_stat		,p_stat_old	,dp_stat_dt	,p_int, T_int, rho_int, rho_old	
	type(field_scalar_cons)	,target	:: div_v	,div_v_int	,ddiv_v_dt	,H	, H_old, R
	type(field_scalar_cons)	,target	:: E_f_int
	type(field_vector_cons)	,target	:: v_int, Y_int, Y_old
	type(field_scalar_flow)	,target	:: F_a, F_b
	type(field_vector_flow)	,target	:: v_f, v_f_old	

	real(dkind)	,dimension(:)	,allocatable	:: concs
	real(dkind)	,dimension(:)	,allocatable	:: Y_rho_int
	!$omp threadprivate(concs,Y_rho_int)

	real(dkind)	,dimension(:)	,allocatable	:: farfield_velocity_array, flame_front_coords
	integer	:: flame_loc_unit, flame_structure_unit

	type subgrid 
		real(dkind), dimension(:,:,:), allocatable :: cells
    end type
    
	type subgrid_vector
		real(dkind), dimension(:,:,:,:), allocatable :: cells
	end type

	type fds_solver
		logical			:: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, hydrodynamics_flag, CFL_condition_flag, all_Neumann_flag
		real(dkind)		:: courant_fraction
		real(dkind)		:: time, time_step, initial_time_step
		integer			:: additional_particles_phases_number, additional_droplets_phases_number
		
		type(viscosity_solver)				:: visc_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(diffusion_solver)				:: diff_solver
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(table_approximated_real_gas)	:: state_eq

!		type(lagrangian_droplets_solver), dimension(:)	    ,allocatable	:: droplets_solver			!# Lagrangian droplets solver
		type(droplets_solver)			, dimension(:)	    ,allocatable	:: droplets_solver			!# Continuum droplets solver
		
		type(lagrangian_particles_solver), dimension(:)	    ,allocatable	:: particles_solver			!# Lagrangian particles solver
!		type(particles_solver)			, dimension(:)	    ,allocatable	:: particles_solver			!# Continuum particles solver
		
		type(computational_domain)					:: domain
		type(thermophysical_properties_pointer)		:: thermo		
		type(chemical_properties_pointer)			:: chem
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		type(field_scalar_cons_pointer)	:: rho		, rho_int		, rho_old		, T				, T_int			, p				, p_int			, v_s			, mol_mix_conc
		type(field_scalar_cons_pointer)	:: E_f		, E_f_prod_chem	, E_f_prod_heat	, E_f_prod_gd	, E_f_prod_visc	, E_f_prod_diff	, E_f_int		, h_s			, gamma
		type(field_scalar_cons_pointer)	:: p_stat	, p_stat_old	, dp_stat_dt	, p_dyn			, div_v			, div_v_int		, ddiv_v_dt		, H				, H_old			, R
		type(field_scalar_cons_pointer)	:: nu		, kappa
		type(field_scalar_flow_pointer)	:: F_a		, F_b
		
		type(field_vector_cons_pointer)	:: v		, v_prod_gd		, v_prod_visc	, v_prod_source	, v_int	
		type(field_vector_cons_pointer)	:: Y		, Y_prod_diff	, Y_prod_chem	, Y_int			, Y_old
		type(field_vector_cons_pointer)	:: D
		type(field_vector_flow_pointer)	:: v_f		, v_f_old
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	::  rho_prod_droplets, E_f_prod_droplets, E_f_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  Y_prod_droplets		
!		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Lagrangian droplets solver
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Continuum droplets solver
		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Lagrangian particles solver
!		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Continuum particles solver		


		real(dkind)	,dimension(:,:,:,:)	,allocatable	:: vorticity, grad_F_a, grad_F_b
		real(dkind)	,dimension(:,:,:)	,allocatable	:: p_old
		
		integer												:: number_of_meshes
		integer			,dimension(:,:,:)	,allocatable	:: subgrid_cons_inner_loop, subgrid_cons_utter_loop
		type(subgrid)	,dimension(:)		,allocatable	:: sub_E, sub_F, sub_R, sub_E_old, sub_bc
        type(subgrid_vector)	,dimension(:)		,allocatable	:: sub_mesh
        real(dkind)				,dimension(:)		,allocatable	:: sub_cells_number
        real(dkind)				,dimension(:)		,allocatable	:: sub_cell_size
		
		real(dkind)				:: rho_0
	contains
		procedure	,private	:: calculate_interm_Y_predictor
		procedure	,private	:: calculate_divergence_v
		procedure	,private	:: calculate_pressure_poisson
		procedure	,private	:: calculate_dynamic_pressure
		procedure	,private	:: calculate_velocity
		procedure	,private	:: calculate_interm_Y_corrector
		procedure	,private	:: CHARM_flux_limiter
		procedure	,private	:: apply_boundary_conditions
		procedure	,private	:: apply_poisson_boundary_conditions
		procedure	,private	:: farfield_values_modifier
		procedure	,private	:: V_cycle
		procedure	,private	:: perturb_velocity_field
		procedure	,private	:: stabilizing_inlet
		procedure	,private	:: stabilizing_inlet_1D
		procedure	,private	:: write_data_table
		procedure	,private	:: igniter
        procedure	,private	:: if_stabilized
		procedure				:: solve_problem
		procedure				:: calculate_time_step
		procedure				:: get_time_step
		procedure				:: get_time
	end type

	interface	fds_solver_c
		module procedure	constructor
	end interface

contains
	
	type(fds_solver)	function constructor(manager,problem_data_io, problem_solver_options)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io
		type(solver_options)					,intent(in)		:: problem_solver_options

		real(dkind)						:: calculation_time		
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
		
		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr		
		
		integer	,dimension(3,2)	:: cons_allocation_bounds		
		integer	,dimension(3,2)	:: cons_utter_loop, cons_inner_loop
		integer	,dimension(3,2)	:: flow_inner_loop, loop
		
		integer				:: species_number

		type(liquid_droplets_phase)     :: droplets_params
		type(solid_particles_phase)     :: particles_params
		integer				:: particles_phase_counter, droplets_phase_counter		
		
		character(len=40)	:: var_name
		
		integer				:: bound_number, plus, sign, bound
		character(len=20)	:: boundary_type_name

		real(dkind)	,dimension(3)	:: cell_size, offset
		real(dkind)					:: x, y
		real(dkind)					:: farfield_velocity
		integer	:: dimensions
		integer	:: i,j,k,dim,dim1
		integer	:: ilb,	jlb, irb, jrb, ilt, jlt, irt, jrt   
		integer	:: mesh
		
		integer :: cells_number
		
		constructor%diffusion_flag		= problem_solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= problem_solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= problem_solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= problem_solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= problem_solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= problem_solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= problem_solver_options%get_CFL_condition_flag()
		constructor%sources_flag		= .false.
		
		constructor%additional_droplets_phases_number	= problem_solver_options%get_additional_droplets_phases_number()		
		constructor%additional_particles_phases_number	= problem_solver_options%get_additional_particles_phases_number()
		
		constructor%domain				= manager%domain
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr		
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%boundary%bc_ptr		=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr

		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr	

		call manager%create_scalar_field(p_dyn		,'pressure_dynamic'						,'p_dyn')
		constructor%p_dyn%s_ptr						=> p_dyn	
		call manager%create_scalar_field(p_stat		,'pressure_static'						,'p_stat')
		constructor%p_stat%s_ptr					=> p_stat		
		call manager%create_scalar_field(p_stat_old	,'pressure_static_old'					,'p_stat_old')
		constructor%p_stat_old%s_ptr				=> p_stat_old			
		call manager%create_scalar_field(div_v		,'velocity_divergence'					,'div_v')
		constructor%div_v%s_ptr						=> div_v			
		call manager%create_scalar_field(H			,'stagnation_energy'					,'H')
		constructor%H%s_ptr							=> H	
		call manager%create_scalar_field(dp_stat_dt	,'pressure_static_change'				,'dp_stat_dt')
		constructor%dp_stat_dt%s_ptr				=> dp_stat_dt	
	
		call manager%create_scalar_field(rho_int	,'density_interm'						,'rho_int')
		constructor%rho_int%s_ptr				=> rho_int
		call manager%create_scalar_field(rho_old	,'density_old'							,'rho_old')
		constructor%rho_old%s_ptr				=> rho_old	
		call manager%create_scalar_field(div_v_int	,'velocity_divergence_interm'			,'div_v_int')
		constructor%div_v_int%s_ptr				=> div_v_int
		call manager%create_scalar_field(ddiv_v_dt	,'velocity_divergence_change'			,'ddiv_v_dt')
		constructor%ddiv_v_dt%s_ptr				=> ddiv_v_dt			
		call manager%create_scalar_field(E_f_int	,'full_energy_interm'					,'E_f_int')
		constructor%E_f_int%s_ptr				=> E_f_int
		call manager%create_scalar_field(T_int		,'temperature_interm'					,'T_int')
		constructor%T_int%s_ptr					=> T_int		
		call manager%create_scalar_field(p_int		,'pressure_interm'						,'p_int')
		constructor%p_int%s_ptr					=> p_int	
		call manager%create_scalar_field(H_old		,'stagnation_energy_old'				,'H_old')
		constructor%H_old%s_ptr					=> H_old	
		call manager%create_scalar_field(R			,'residual'								,'R')
		constructor%R%s_ptr							=> R			
		
		call manager%create_vector_field(v_int		,'velocity_interm'						,'v_int'	,'spatial')
		constructor%v_int%v_ptr					=> v_int
		call manager%create_vector_field(Y_int		,'specie_molar_concentration_interm'	,'Y_int'	,'chemical')
		constructor%Y_int%v_ptr					=> Y_int
		call manager%create_vector_field(Y_old		,'specie_molar_concentration_old'		,'Y_old'	,'chemical')
		constructor%Y_old%v_ptr					=> Y_old
		
		call manager%create_scalar_field(F_a	,'F_a'				,'F_a')
		constructor%F_a%s_ptr 	=> F_a
		call manager%create_scalar_field(F_b	,'F_b'				,'F_b')
		constructor%F_b%s_ptr 	=> F_b
		
		call manager%create_vector_field(v_f	,'velocity_flow'						,'v_f'		,'spatial')
		constructor%v_f%v_ptr => v_f
		call manager%create_vector_field(v_f_old,'velocity_flow_old'					,'v_f_old'	,'spatial')
		constructor%v_f_old%v_ptr => v_f_old
		
		constructor%state_eq	=	table_approximated_real_gas_c(manager)
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'sensible_enthalpy')
		constructor%h_s%s_ptr			=> scal_ptr%s_ptr	
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'adiabatic_index')
		constructor%gamma%s_ptr			=> scal_ptr%s_ptr			
		
		!if(constructor%viscosity_flag) then
		!	constructor%visc_solver			= viscosity_solver_c(manager)
		!	call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity')
		!	constructor%E_f_prod_visc%s_ptr			=> scal_ptr%s_ptr
		!	call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
		!	constructor%v_prod_visc%v_ptr			=> vect_ptr%v_ptr
		!	call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
		!	constructor%nu%s_ptr					=> scal_ptr%s_ptr			
		!end if
        
 		if(constructor%viscosity_flag) then
			constructor%visc_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity_FDS')
			constructor%E_f_prod_visc%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_ptr%v_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'viscosity')
			constructor%nu%s_ptr					=> scal_ptr%s_ptr			
		end if

		if (constructor%heat_trans_flag) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr		=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'thermal_conductivity')
			constructor%kappa%s_ptr				=> scal_ptr%s_ptr
		end if

		if (constructor%diffusion_flag) then
			constructor%diff_solver			= diffusion_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_diffusion')
			constructor%E_f_prod_diff%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_ptr%v_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'diffusivity')
			constructor%D%v_ptr						=> vect_ptr%v_ptr
		end if

		if(constructor%additional_droplets_phases_number /= 0) then
			allocate(constructor%droplets_solver(constructor%additional_droplets_phases_number))
			call constructor%droplets_solver(1)%pre_constructor(constructor%additional_droplets_phases_number)
			allocate(constructor%rho_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%E_f_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%v_prod_droplets(constructor%additional_droplets_phases_number))
			allocate(constructor%Y_prod_droplets(constructor%additional_droplets_phases_number))
			do droplets_phase_counter = 1, constructor%additional_droplets_phases_number
				droplets_params = problem_solver_options%get_droplets_params(droplets_phase_counter)
!				constructor%droplets_solver(droplets_phase_counter)	= lagrangian_droplets_solver_c(manager, droplets_params, droplets_phase_counter)		!# Lagrangian droplets solver
				constructor%droplets_solver(droplets_phase_counter)	= droplets_solver_c(manager, droplets_params, droplets_phase_counter)					!# Continuum droplets solver
				write(var_name,'(A,I2.2)') 'energy_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_droplets', droplets_phase_counter						
!				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)								!# Lagrangian droplets solver
!				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr						!# Lagrangian droplets solver
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum droplets solver								
				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum droplets solver
				write(var_name,'(A,I2.2)') 'concentration_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%Y_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr                
			end do		
		end if   
		
		if(constructor%additional_particles_phases_number /= 0) then
			allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
			call constructor%particles_solver(1)%pre_constructor(constructor%additional_particles_phases_number)
			allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
			do particles_phase_counter = 1, constructor%additional_particles_phases_number
				particles_params = problem_solver_options%get_particles_params(particles_phase_counter)
				constructor%particles_solver(particles_phase_counter)	= lagrangian_particles_solver_c(manager, particles_params, particles_phase_counter)	!# Lagrangian particles solver
!				constructor%particles_solver(particles_phase_counter)	= particles_solver_c(manager, particles_params, particles_phase_counter)			!# Continuum particles solver
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter						
				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)								!# Lagrangian particles solver
				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr						!# Lagrangian particles solver
!				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum particles solver								
!				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum particles solver
			end do		
		end if
		
		if (constructor%reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_ptr%v_ptr
		end if

		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()		
		
		allocate(constructor%vorticity(		3						, &
											cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))		

		allocate(constructor%p_old(			cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))													

		allocate(constructor%grad_F_a(		3						, &
											cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
				
		allocate(constructor%grad_F_b(		3						, &
											cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
											cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
											cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		constructor%vorticity= 0.0_dkind
		constructor%p_old	 = 0.0_dkind
		constructor%grad_F_a = 0.0_dkind									
		constructor%grad_F_b = 0.0_dkind

		cons_utter_loop	= manager%domain%get_local_utter_cells_bounds()	
		cons_inner_loop = manager%domain%get_local_inner_cells_bounds()	
		flow_inner_loop	= manager%domain%get_local_inner_faces_bounds()	
		
		problem_data_io				= data_io_c(manager,calculation_time)									
		
		call problem_data_io%input_all_data()
		
		do dim = 1, dimensions		

			loop = flow_inner_loop

			do dim1 = 1,dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2) 
			do i = loop(1,1),loop(1,2)
			
				do dim1 = 1, dimensions
					constructor%v_f%v_ptr%pr(dim1)%cells(dim,i,j,k) = 0.5_dkind * (constructor%v%v_ptr%pr(dim1)%cells(i,j,k) + constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) )
				end do			
			
			end do
			end do
			end do
				
		end do
		
		if(problem_data_io%get_load_counter() == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
			call problem_data_io%add_io_scalar_cons_field(constructor%div_v)
		
			constructor%p_dyn%s_ptr%cells		= 0.0_dkind      
			constructor%p_stat%s_ptr%cells		= constructor%p%s_ptr%cells	
			constructor%v_f%v_ptr%pr(1)%cells	= 0.0_dkind 
			
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
			
			if(constructor%additional_droplets_phases_number /= 0) then
				do droplets_phase_counter = 1, constructor%additional_droplets_phases_number
					call constructor%droplets_solver(droplets_phase_counter)%set_initial_distributions()
				end do
			end if  
			
			if(constructor%additional_particles_phases_number /= 0) then
				do particles_phase_counter = 1, constructor%additional_particles_phases_number
					call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
				end do
			end if 			
			
		end if		
		dimensions		= manager%domain%get_domain_dimensions()

		cell_size						= constructor%mesh%mesh_ptr%get_cell_edges_length()
		
		constructor%time				= calculation_time
		constructor%initial_time_step	= problem_solver_options%get_initial_time_step()
		constructor%time_step			= constructor%initial_time_step

		constructor%rho_0				= constructor%rho%s_ptr%cells(1,1,1)!constructor%rho%s_ptr%cells(1,1,1)
		print *,constructor%rho_0
	
		species_number = manager%chemistry%chem_ptr%species_number
		
		do bound = 1, size(constructor%boundary%bc_ptr%boundary_types)
			boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound)%get_type_name()
			if (boundary_type_name == 'inlet') then
				farfield_velocity = constructor%boundary%bc_ptr%boundary_types(bound)%get_farfield_velocity()
			end if
		end do		

		allocate(farfield_velocity_array(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))
		allocate(flame_front_coords(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))
		
		flame_front_coords = 0.0
		farfield_velocity_array = farfield_velocity

		!$omp parallel
		allocate(concs(species_number))
		allocate(Y_rho_int(species_number))
		concs		= 0.0_dkind
		Y_rho_int	= 0.0_dkind
		!$omp end parallel
		
		cells_number = cons_inner_loop(1,2) - cons_inner_loop(1,1) + 1
		do dim = 2, dimensions
			cells_number = min(cells_number, cons_inner_loop(dim,2) - cons_inner_loop(dim,1) + 1)
		end do
			
		constructor%number_of_meshes = 1 
		
		do while(mod(cells_number,2) == 0)
			cells_number = cells_number /2
			constructor%number_of_meshes = constructor%number_of_meshes + 1
		end do
		
		allocate(constructor%sub_E(0:constructor%number_of_meshes-1),constructor%sub_E_old(0:constructor%number_of_meshes-1),constructor%sub_F(0:constructor%number_of_meshes-1),constructor%sub_R(0:constructor%number_of_meshes-1),constructor%sub_bc(0:constructor%number_of_meshes-1),constructor%sub_mesh(0:constructor%number_of_meshes-1))

		allocate(constructor%subgrid_cons_utter_loop(0:constructor%number_of_meshes-1,3,2), constructor%subgrid_cons_inner_loop(0:constructor%number_of_meshes-1,3,2))
		
        allocate(constructor%sub_cells_number(0:constructor%number_of_meshes-1))
        allocate(constructor%sub_cell_size(0:constructor%number_of_meshes-1))
        
		constructor%subgrid_cons_inner_loop = 1
		constructor%subgrid_cons_utter_loop = 1
        
		do mesh = 0,constructor%number_of_meshes-1

            cells_number = 1
            
			do dim = 1, dimensions
				constructor%subgrid_cons_utter_loop(mesh,dim,1) = cons_utter_loop(dim,1)
				constructor%subgrid_cons_utter_loop(mesh,dim,2) = (cons_utter_loop(dim,2)-1)/2**mesh + 1
				constructor%subgrid_cons_inner_loop(mesh,dim,1) = constructor%subgrid_cons_utter_loop(mesh,dim,1) + 1
				constructor%subgrid_cons_inner_loop(mesh,dim,2) = constructor%subgrid_cons_utter_loop(mesh,dim,2) - 1
			end do

			allocate(constructor%sub_E(mesh)%cells(	constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))
													
			allocate(constructor%sub_F(mesh)%cells(	constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))
													
			allocate(constructor%sub_R(mesh)%cells(	constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
													constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))
													
			allocate(constructor%sub_E_old(mesh)%cells(	constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))	
														
			allocate(constructor%sub_bc(mesh)%cells(	constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))
            
			allocate(constructor%sub_mesh(mesh)%cells(	dimensions, & 
														constructor%subgrid_cons_utter_loop(mesh,1,1):constructor%subgrid_cons_utter_loop(mesh,1,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,2,1):constructor%subgrid_cons_utter_loop(mesh,2,2)	, &
														constructor%subgrid_cons_utter_loop(mesh,3,1):constructor%subgrid_cons_utter_loop(mesh,3,2)))            
			
														

			constructor%sub_E(mesh)%cells		= 0.0_dkind
			constructor%sub_F(mesh)%cells		= 0.0_dkind
			constructor%sub_R(mesh)%cells		= 0.0_dkind
			constructor%sub_E_old(mesh)%cells	= 0.0_dkind
			constructor%sub_bc(mesh)%cells		= 0.0_dkind
            constructor%sub_mesh(mesh)%cells	= 0.0_dkind
            
            do dim = 1, dimensions
				cells_number = cells_number * (constructor%subgrid_cons_inner_loop(mesh,dim,2) - constructor%subgrid_cons_inner_loop(mesh,dim,1) + 1)
            end do
		
            constructor%sub_cells_number(mesh) = cells_number
            constructor%sub_cell_size(mesh) = cell_size(1)*2**mesh
		end do	
		
		constructor%sub_bc(0)%cells = constructor%boundary%bc_ptr%bc_markers
		
		constructor%all_Neumann_flag = .true.
		
		do mesh = 0,constructor%number_of_meshes-2

			do k = constructor%subgrid_cons_utter_loop(mesh,3,1),constructor%subgrid_cons_utter_loop(mesh,3,2)
			do j = constructor%subgrid_cons_utter_loop(mesh,2,1),constructor%subgrid_cons_utter_loop(mesh,2,2)
			do i = constructor%subgrid_cons_utter_loop(mesh,1,1),constructor%subgrid_cons_utter_loop(mesh,1,2)
				if(constructor%sub_bc(mesh)%cells(i,j,k) /= 0) then

					if (dimensions == 1) then
				
						ilb = int(i/2)
						irb = int(i/2)+1
					
						if (mod(i,2) == 0) then
							constructor%sub_bc(mesh+1)%cells(ilb,j,k) = constructor%sub_bc(mesh)%cells(i,j,k)  
						end if			
						if (mod(i,2) /= 0) then
							constructor%sub_bc(mesh+1)%cells(irb,j,k) = constructor%sub_bc(mesh)%cells(i,j,k)
                        end if	
                         
					end if					
				
					if (dimensions == 2) then

						ilb = int(i/2)
						jlb = int(j/2)
							
						irb = int(i/2)+1
						jrb = int(j/2)
							
						ilt = int(i/2)
						jlt = int(j/2)+1
							
						irt = int(i/2)+1
						jrt = int(j/2)+1
						
						if ((mod(i,2) == 0).and.(mod(j,2) /= 0)) then
							if (constructor%sub_bc(mesh)%cells(i,j,k) > constructor%sub_bc(mesh+1)%cells(ilt,jlt,k)) then 
								constructor%sub_bc(mesh+1)%cells(ilt,jlt,k) = constructor%sub_bc(mesh)%cells(i,j,k)
							end if
						end if
						if ((mod(i,2) /= 0).and.(mod(j,2) == 0)) then
							if (constructor%sub_bc(mesh)%cells(i,j,k) > constructor%sub_bc(mesh+1)%cells(irb,jrb,k)) then
								constructor%sub_bc(mesh+1)%cells(irb,jrb,k) = constructor%sub_bc(mesh)%cells(i,j,k)
							end if
						end if
						if ((mod(i,2) == 0).and.(mod(j,2) == 0)) then
							if (constructor%sub_bc(mesh)%cells(i,j,k) > constructor%sub_bc(mesh+1)%cells(ilb,jlb,k)) then
								constructor%sub_bc(mesh+1)%cells(ilb,jlb,k) = constructor%sub_bc(mesh)%cells(i,j,k)
							end if
						end if			
						if ((mod(i,2) /= 0).and.(mod(j,2) /= 0)) then
							if (constructor%sub_bc(mesh)%cells(i,j,k) > constructor%sub_bc(mesh+1)%cells(irt,jrt,k)) then
								constructor%sub_bc(mesh+1)%cells(irt,jrt,k) = constructor%sub_bc(mesh)%cells(i,j,k)
							end if
						end if	

					end if
					
					bound_number = constructor%boundary%bc_ptr%bc_markers(i,j,k)
					
					if (mesh == 0) then
						boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound_number)%get_type_name()
									
						if (boundary_type_name == 'outlet') then
							constructor%all_Neumann_flag = .false.
						end if
					end if
							
				end if
			end do
			end do
			end do


            constructor%sub_mesh(0)%cells =  constructor%mesh%mesh_ptr%mesh
            
			do k = constructor%subgrid_cons_inner_loop(mesh+1,3,1),constructor%subgrid_cons_inner_loop(mesh+1,3,2)
			do j = constructor%subgrid_cons_inner_loop(mesh+1,2,1),constructor%subgrid_cons_inner_loop(mesh+1,2,2)
			do i = constructor%subgrid_cons_inner_loop(mesh+1,1,1),constructor%subgrid_cons_inner_loop(mesh+1,1,2)            
				offset = 1
				do dim = 1, dimensions
					offset(dim) = constructor%subgrid_cons_inner_loop(mesh,dim,1) + 2 * ((i-1) * I_m(dim,1) + (j-1) * I_m(dim,2) + (k-1) * I_m(dim,3))
				end do
				
				do dim = 1, dimensions
					constructor%sub_mesh(mesh+1)%cells(dim,i,j,k) = 0.5_dkind*(constructor%sub_mesh(mesh)%cells(dim,offset(1),offset(2),offset(3)) + constructor%sub_mesh(mesh)%cells(dim,offset(1)+1,offset(2),offset(3)))
				end do	
            
			end do
			end do
            end do    
            
            do dim = 1, dimensions
                i = constructor%subgrid_cons_utter_loop(mesh+1,1,1) * I_m(dim,1) + (1-I_m(dim,1))
                j = constructor%subgrid_cons_utter_loop(mesh+1,2,1) * I_m(dim,2) + (1-I_m(dim,2))
                k = constructor%subgrid_cons_utter_loop(mesh+1,3,1) * I_m(dim,3) + (1-I_m(dim,3))
                
				constructor%sub_mesh(mesh+1)%cells(dim,i,j,k) = constructor%sub_mesh(mesh+1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - cell_size(1)*2**(mesh+1)
                
                i = constructor%subgrid_cons_utter_loop(mesh+1,1,2) * I_m(dim,1) + (1-I_m(dim,1))
                j = constructor%subgrid_cons_utter_loop(mesh+1,2,2) * I_m(dim,2) + (1-I_m(dim,2))
                k = constructor%subgrid_cons_utter_loop(mesh+1,3,2) * I_m(dim,3) + (1-I_m(dim,3))
                
                constructor%sub_mesh(mesh+1)%cells(dim,i,j,k) = constructor%sub_mesh(mesh+1)%cells(dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + cell_size(1)*2**(mesh+1)
            end do
            
        end do
		
		!!!!!######
		!constructor%sub_bc(5)%cells(2,0,1) = 3.0
		
		open(newunit = flame_loc_unit, file = 'av_flame_data.dat', status = 'replace', form = 'formatted')

		continue
		
	end function

	subroutine solve_problem(this,iteration,stop_flag)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: iteration
        logical				,intent(inout)	:: stop_flag
		
		integer	:: droplets_phase_counter, particles_phase_counter
		integer	:: specie

		logical	:: perturbed_velocity_field, stabilizing_inlet, ignite, stabilized
		
		perturbed_velocity_field	= .false.
		stabilizing_inlet			= .false.
		ignite						= .false.
		
		this%time = this%time + this%time_step		

		if (ignite) then
			call this%igniter(this%time)
		end if
		
!		call this%farfield_values_modifier(this%time)
		
!		if (iteration > 10) perturbed_velocity_field = .false.
		if (perturbed_velocity_field)	call this%perturb_velocity_field(this%time_step)

		if (this%reactive_flag)		call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
		if (this%viscosity_flag)	call this%visc_solver%solve_viscosity(this%time_step)
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)				
		
		if(this%additional_droplets_phases_number /= 0) then
			do droplets_phase_counter = 1, this%additional_droplets_phases_number
!				call this%droplets_solver(droplets_phase_counter)%droplets_solve(this%time_step)				!# Lagrangian droplets solver
				call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_main(this%time)		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_euler_step_v_E(this%time_step)		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_lagrange_step(this%time_step)		!# Continuum droplets solver
				call this%droplets_solver(droplets_phase_counter)%droplets_final_step(this%time_step)			!# Continuum droplets solver		
			end do		
		end if  	
		
		if(this%additional_particles_phases_number /= 0) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
				call this%particles_solver(particles_phase_counter)%particles_solve(this%time_step)				!# Lagrangian particles solver
				!call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_main()			!# Continuum particles solver
				!call this%particles_solver(particles_phase_counter)%particles_euler_step_v_E(this%time_step)	!# Continuum particles solver
				!call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum particles solver
				!call this%particles_solver(particles_phase_counter)%particles_lagrange_step(this%time_step)	!# Continuum particles solver
				!call this%particles_solver(particles_phase_counter)%particles_final_step(this%time_step)		!# Continuum particles solver		
			end do		
		end if
		
		call this%calculate_interm_Y_predictor(this%time_step)
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.true.)		
		call this%apply_boundary_conditions(this%time_step,predictor=.true.)		
		call this%calculate_divergence_v		(this%time_step,predictor=.true.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.true.)
		call this%calculate_velocity			(this%time_step,predictor=.true.)
		
		if (this%viscosity_flag)	call this%visc_solver%solve_viscosity(this%time_step)

		call this%calculate_interm_Y_corrector(this%time_step)
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.false.)
		call this%apply_boundary_conditions(this%time_step,predictor=.false.)

		call this%calculate_divergence_v		(this%time_step,predictor=.false.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.false.)
		call this%calculate_velocity			(this%time_step,predictor=.false.)		

		call this%calculate_time_step()
		
!		call this%chem_kin_solver%write_chemical_kinetics_table('15_pcnt_H2-Air_table(T).dat')

!        call this%if_stabilized(this%time, stabilized)
!        if (stabilized) stop
        
		if (stabilizing_inlet) then
			call this%stabilizing_inlet_1D(this%time, stabilized)
			if (stabilized) then
               stop_flag = .true.
!				call this%chem_kin_solver%write_chemical_kinetics_table('29.5_pcnt_H2-Air_table(T).dat')
!				call this%write_data_table('H2-Air_initials.dat')
				stop
            end if
		end if
         
		!call this%state_eq%check_conservation_laws()

	end subroutine

	subroutine calculate_interm_Y_predictor(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		
		real(dkind)	:: B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dkind)	:: spec_summ
		
		real(dkind)	,dimension(3)	:: cell_size		
		real(dkind)	:: rhs
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		integer	:: bound_number
		integer :: i,j,k,dim,spec,droplets_phase_counter
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		coordinate_system	= this%domain%get_coordinate_system_name()

		associate (	rho				=> this%rho%s_ptr			, &
					rho_int			=> this%rho_int%s_ptr		, &
					v_f				=> this%v_f%v_ptr			, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_droplets	=> this%Y_prod_droplets	, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)
					
			!$omp parallel default(none)  private(i,j,k,dim,spec,spec_summ,rhs,flux_right,flux_left,lame_coeffs) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(v_f,Y_prod_diff,Y_prod_chem,Y,Y_old,Y_int,rho,rho_old,rho_int,bc,cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,mesh,time_step)					
			!$omp do collapse(3) schedule(guided)	
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	

					lame_coeffs		= 1.0_dkind				
				
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dkind
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
					end select					
				
					do spec = 1, species_number
					
						Y_rho_int(spec) = rho%cells(i,j,k)*Y%pr(spec)%cells(i,j,k) 
						
						rhs = 0.0_dkind
						
						do dim = 1, dimensions
													
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
									flux_right = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								else
									flux_right = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end if
							end if
						
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
								flux_left = this%CHARM_flux_limiter(rho%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3))*		&
																	Y%pr(spec)%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)),	&
																	v_f%pr(dim)%cells(dim,i,j,k))
							else
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dkind) then
									flux_left = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								else
									flux_left = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								end if
							end if

							rhs = rhs -(	flux_right * v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) &
										-	flux_left  * v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2)
						
							continue			
						end do
						
						Y_rho_int(spec) = Y_rho_int(spec)  + rhs* time_step
						
						if (this%diffusion_flag)	Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_diff%pr(spec)%cells(i,j,k) * time_step	![kg/m^3]
						if (this%reactive_flag)		Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_chem%pr(spec)%cells(i,j,k) * time_step	![kg/m^3]		
						
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
								Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_droplets(droplets_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step		![kg/m^3]
							end do		
						end if						
						
					end do	
					
					rho_int%cells(i,j,k)	= sum(Y_rho_int)
					
					do spec = 1, species_number
						Y_int%pr(spec)%cells(i,j,k) = Y_rho_int(spec) / rho_int%cells(i,j,k) 
					end do
				end if
			end do
			end do
			end do				
			!$omp end do
			
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					spec_summ = 0.0_dkind
					do spec = 1,species_number
						spec_summ = spec_summ + max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind)
					end do
					do spec = 1,species_number
						Y_int%pr(spec)%cells(i,j,k) = max(Y_int%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ
					end do
				end if
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then						
					rho_old%cells(i,j,k)	= rho%cells(i,j,k)
					rho%cells(i,j,k)		= rho_int%cells(i,j,k)			
					do spec = 1, species_number
						Y_old%pr(spec)%cells(i,j,k)	= Y%pr(spec)%cells(i,j,k)
						Y%pr(spec)%cells(i,j,k)		= Y_int%pr(spec)%cells(i,j,k)			
					end do						
				end if
			end do
			end do
			end do
			!$omp end do
			!$omp end parallel
			
		!	print *, '1',Y_int%pr(1)%cells(50,20,1),rho_int%cells(50,20,1)
			
			continue
			
		end associate

	end subroutine
	
	subroutine calculate_divergence_v(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
	
		real(dkind)	:: flux_left, flux_right, specie_enthalpy, mixture_cp, dp_dt, D_sum, P_sum, U_sum, div_sum, average_molar_mass, mol_mix_conc
		
		real(dkind)					:: cell_volume
		real(dkind)	,dimension(3)	:: cell_size, cell_surface_area
			
		real(dkind)					:: energy_source = 1.0e06_dkind
		real(dkind)	,save			:: time
		integer		,save			:: iter = 1

		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system	
		
		integer	:: dimensions, species_number, cells_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number, plus, sign
		integer :: i,j,k,dim,spec,droplets_phase_counter
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		cell_volume		= this%mesh%mesh_ptr%get_cell_volume()
		
		
		associate (	div_v_int		=> this%div_v_int%s_ptr		, &
					T				=> this%T%s_ptr				, &
					rho				=> this%rho%s_ptr			, &
					p_stat			=> this%p_stat%s_ptr		, &
					gamma			=> this%gamma%s_ptr			, &
					h_s				=> this%h_s%s_ptr			, &
					dp_stat_dt		=> this%dp_stat_dt%s_ptr	, &
					v_f				=> this%v_f%v_ptr			, &
					Y				=> this%Y%v_ptr				, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &	
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &	
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_gd 	=> this%E_f_prod_gd%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff	=> this%E_f_prod_diff%s_ptr	, &	
					E_f_prod_droplets	=> this%E_f_prod_droplets, &
					thermo			=> this%thermo%thermo_ptr	, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

			D_sum = 0.0_dkind
			P_sum = 0.0_dkind		
			U_sum = 0.0_dkind
			
			div_sum = 0.0_dkind
			
			time = time + 0.5_dkind*time_step
			
			if(time <= 1.0e-02_dkind) then
				energy_source = 2.0e02_dkind
			else
				energy_source = 0.0_dkind
			end if
			
		!	print *, '2', div_v_int%cells(50,20,1),T%cells(50,20,1),D_sum,P_sum
			
			!$omp parallel default(none)  private(flux_right,flux_left,i,j,k,dim,spec,mixture_cp,specie_enthalpy,plus,sign,bound_number,cell_volume,cell_surface_area,lame_coeffs, average_molar_mass, mol_mix_conc) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(D_sum,P_sum,U_sum,div_v_int,p_stat,dp_stat_dt,v_f,h_s,E_f_prod_visc,E_f_prod_heat,E_f_prod_diff,E_f_prod_chem,Y_prod_diff,Y_prod_chem,T,Y,rho,energy_source,bc,thermo,cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,mesh)					
			!$omp do collapse(3) schedule(guided)	reduction(+:D_sum,P_sum)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
				
					cell_volume	= mesh%get_cell_volume()
					lame_coeffs	= 1.0_dkind

					select case(coordinate_system)
						case ('cartesian')	

						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
							cell_volume			= cell_volume * mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
							cell_volume			= cell_volume * mesh%mesh(1,i,j,k)**2
					end select					

					div_v_int%cells(i,j,k) = 0.0_dkind
					
					!if ((i > cons_inner_loop(1,1) + 1).and.(j > cons_inner_loop(2,1) + 1).and.(i < cons_inner_loop(1,2) - 1).and.(j < cons_inner_loop(2,2) - 1)) then
					if (this%viscosity_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_visc%cells(i,j,k)	![J/m^3/s]
					if (this%heat_trans_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_heat%cells(i,j,k)	![J/m^3/s]
					if (this%diffusion_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_diff%cells(i,j,k)	![J/m^3/s]
					if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_chem%cells(i,j,k)	![J/m^3/s]
					!end if

					average_molar_mass = 0.0_dkind
					do spec = 1,species_number
						average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
					end do
				
					mol_mix_conc		= 1.0_dkind / average_molar_mass

					concs = 0.0_dkind
					do spec = 1,species_number
						concs(spec)				= Y%pr(spec)%cells(i,j,k) *  mol_mix_conc / thermo%molar_masses(spec)
					end do		

					if (this%additional_droplets_phases_number /= 0) then
						do droplets_phase_counter = 1, this%additional_droplets_phases_number
							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + E_f_prod_droplets(droplets_phase_counter)%s_ptr%cells(i,j,k)	![J/m^3/s]
						end do		
					end if
					
					do dim = 1,dimensions
					
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
							flux_right = this%CHARM_flux_limiter(	rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	h_s%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),		&
																	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
						else
							if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
								flux_right = rho%cells(i,j,k) * h_s%cells(i,j,k)
							else
								flux_right = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * h_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
							end if
						end if

						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
							flux_left = this%CHARM_flux_limiter(	rho%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3))*		&
																	h_s%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)),		&
																	v_f%pr(dim)%cells(dim,i,j,k))
							continue
						else
							if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dkind) then
								flux_left = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * h_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
							else
								flux_left = rho%cells(i,j,k) * h_s%cells(i,j,k)
							end if
						end if					
						
						div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k)  -  (	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* lame_coeffs(dim,3) * (flux_right	-  rho%cells(i,j,k) * h_s%cells(i,j,k)) / cell_size(1)	&
																			  -	v_f%pr(dim)%cells(dim,i,j,k)									* lame_coeffs(dim,1) * (flux_left	-  rho%cells(i,j,k) * h_s%cells(i,j,k)) / cell_size(1)) / lame_coeffs(dim,2)	
						continue
					end do	
				
					mixture_cp				= thermo%calculate_mixture_cp(T%cells(i,j,k), concs)
					
					div_v_int%cells(i,j,k)	= div_v_int%cells(i,j,k) / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc)
					
					do spec = 1, species_number
					
						specie_enthalpy = (thermo%calculate_specie_enthalpy(T%cells(i,j,k),spec) - thermo%calculate_specie_enthalpy(298.15_dkind,spec))  / thermo%molar_masses(spec)
					
						do dim = 1, dimensions
					
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
									flux_right = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								else
									flux_right = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end if
							end if	
						
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
								flux_left = this%CHARM_flux_limiter(rho%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3))*		&
																	Y%pr(spec)%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)),	&
																	v_f%pr(dim)%cells(dim,i,j,k))
							else
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dkind) then
									flux_left = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								else
									flux_left = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								end if
							end if					

							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) *  mixture_cp * T%cells(i,j,k) / mol_mix_conc)) * &
																			  (	-  (	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* lame_coeffs(dim,3) * (flux_right	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)	&
																					 -	v_f%pr(dim)%cells(dim,i,j,k)									* lame_coeffs(dim,1) * (flux_left	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)) / lame_coeffs(dim,2))
							continue													 
						end do							
						
						
						if (this%diffusion_flag) 	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc)) * &
																		  (	Y_prod_diff%pr(spec)%cells(i,j,k))

						if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc)) * &
																		  (	Y_prod_chem%pr(spec)%cells(i,j,k))
					end do
					
					D_sum = D_sum + div_v_int%cells(i,j,k) * cell_volume
					P_sum = P_sum + (1.0_dkind / p_stat%cells(i,j,k) - 1.0_dkind / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc)) * cell_volume 
				end if
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(guided)	reduction(+:U_sum)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
					do dim = 1,dimensions															 
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								cell_surface_area	= mesh%get_cell_surface_area()
								select case(coordinate_system)
									case ('cartesian')	
										cell_surface_area	= cell_surface_area
									case ('cylindrical')
										! x -> r, y -> z
										if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * mesh%mesh(1,i,j,k)									
										if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * mesh%mesh(1,i,j,k)			! - 0.5_dkind*cell_size(1)
									case ('spherical')
										! x -> r
										if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))**2		! - 0.5_dkind*cell_size(1)
								end select	
								U_sum = U_sum + sign * v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) * cell_surface_area(dim)
							end if
						end do	
					end do
				end if
			end do
			end do
			end do
			!$omp end do			

			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)			
				if(bc%bc_markers(i,j,k) == 0) then	
					
					average_molar_mass = 0.0_dkind
					do spec = 1,species_number
						average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
                    end do

					mol_mix_conc	= 1.0_dkind / average_molar_mass				
				
                    concs = 0.0_dkind
					do spec = 1,species_number
						concs(spec)				= Y%pr(spec)%cells(i,j,k) *  mol_mix_conc / thermo%molar_masses(spec)
					end do	    
      
					mixture_cp		= thermo%calculate_mixture_cp(T%cells(i,j,k), concs)
				
					dp_stat_dt%cells(i,j,k) = (D_sum - U_sum)/ P_sum

					if (this%all_Neumann_flag) then	
						div_v_int%cells(i,j,k)  = div_v_int%cells(i,j,k) - (1.0_dkind / p_stat%cells(i,j,k) - 1.0_dkind / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc))* dp_stat_dt%cells(i,j,k)
					end if

				end if
			end do
			end do
			end do			
			!$omp end do
			!$omp end parallel

		
			iter = iter + 1
		!	print *, '2', div_v_int%cells(50,20,1),T%cells(50,20,1),D_sum,P_sum
			
		continue
		
		end associate
		
		
	end subroutine	
	
	subroutine calculate_pressure_poisson(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		
		real(dkind)	:: H_center, H_left, H_right, F_a_left, F_a_right, F_b_left, F_b_right, sub_F_summ, r_summ
		real(dkind)	:: H_residual, H_max, H_max_old, H_average, residual, a_norm_init, a_norm, a_norm_prev, a_normfinal, H_summ, sum_ddiv_v_dt, grad_F_a_summ, grad_F_b_summ
		real(dkind)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		integer		:: r_i, r_j, r_k

		real(dkind)	,dimension(3)	:: cell_size		
		real(dkind)	:: time_step_adj, beta, beta_old, spectral_radii_jacobi, edges_number
		real(dkind), save	:: time = 0.0_dkind
		
		integer	:: dimensions, iterations, iterations_number
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop, subgrid_cons_inner_loop
		integer	,dimension(3,2)	:: loop
		integer	,dimension(3)	:: offset
		
		character(len=20)		:: boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3
		real(dkind)				:: bc_coeff
		
		integer	:: droplets_phase_counter, particles_phase_counter
		integer	:: sign, bound_number, bound_number1, bound_number2, bound_number3
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		integer			:: poisson_iteration, pressure_iteration, overall_poisson_iteration, mesh_iter, v_cycle_iteration
		integer			:: nu_0, nu_1, nu_2, V_cycle_depth
		real(dkind)		:: tolerance
		
		logical			:: converged = .false., v_cycle_converged = .false.
		logical			:: pressure_converged = .false.
		
		dimensions		= this%domain%get_domain_dimensions()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		subgrid_cons_inner_loop	= this%subgrid_cons_inner_loop(1,:,:)
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		time_step_adj	=  0.2_dkind*cell_size(1)**2								!	0.025_dkind*(3.0_dkind * cell_size(1)**2) !0.1_dkind*cell_size(1)**2
		
		if (predictor) time = time + time_step

		associate (	p				=> this%p%s_ptr				, &
					p_stat			=> this%p_stat%s_ptr		, &
					p_dyn			=> this%p_dyn%s_ptr			, &
					rho_old			=> this%rho_old%s_ptr		, &
					rho_int			=> this%rho_int%s_ptr		, &
					div_v_int		=> this%div_v_int%s_ptr		, &
					ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					H				=> this%H%s_ptr				, &
					H_old			=> this%H_old%s_ptr			, &
					R				=> this%R%s_ptr				, &
					vorticity		=> this%vorticity			, &
					grad_F_a		=> this%grad_F_a			, &
					grad_F_b		=> this%grad_F_b			, &
					v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					v_prod_droplets	=> this%v_prod_droplets		, &
					v_prod_particles=> this%v_prod_particles	, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr		, &
					
					sub_F			=> this%sub_F				, &
					sub_R			=> this%sub_R				, &
					sub_E			=> this%sub_E				, &
					sub_E_old		=> this%sub_E_old)

		
			ddiv_v_dt%cells	= 0.0_dkind	
			sum_ddiv_v_dt	= 0.0_dkind
			vorticity		= 0.0_dkind	
			F_a%cells		= 0.0_dkind		
			F_b%cells		= 0.0_dkind		

			grad_F_a_summ	= 0.0_dkind

			
			!$omp parallel default(none)  private(i,j,k,dim,dim1,dim2,loop,lame_coeffs,farfield_velocity,sign,bound_number,bound_number1,bound_number2,bound_number3,plus,boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3,bc_coeff) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(ddiv_v_dt,sum_ddiv_v_dt,div_v_int,vorticity,v_f,v_f_old,F_a,grad_F_a, grad_F_a_summ, rho_old,rho_int,v_prod_visc,bc,cons_inner_loop,cons_utter_loop,flow_inner_loop,dimensions,predictor,cell_size,coordinate_system,mesh,time_step)					
			!$omp do collapse(3) schedule(guided)	reduction(+:sum_ddiv_v_dt) 
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
				
					lame_coeffs		= 1.0_dkind				
		
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dkind
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
					end select					
				
					if (predictor) then
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / time_step
						
						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k) - (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3)  - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2) / time_step
						end do
					else
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / (0.5_dkind * time_step)
						
						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k)  - 0.5_dkind * ( (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)	&
																							+ (v_f_old%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f_old%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)) / (0.5_dkind * time_step)
						end do
					end if

					sum_ddiv_v_dt = sum_ddiv_v_dt + ddiv_v_dt%cells(i,j,k)
				end if
			end do
			end do
			end do	
			!$omp end do

			!if (this%all_Neumann_flag) then
			!!$omp do collapse(3) schedule(guided)
			!do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			!do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			!do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			!	if(bc%bc_markers(i,j,k) == 0) then
			!		ddiv_v_dt%cells(i,j,k)	=	ddiv_v_dt%cells(i,j,k) - sum_ddiv_v_dt/this%sub_cells_number(0)
			!	end if
			!end do
			!end do
			!end do		
			!!$omp end do			
			!end if
			
			!$omp do collapse(3) schedule(guided)	
			do k = flow_inner_loop(3,1),flow_inner_loop(3,2)
			do j = flow_inner_loop(2,1),flow_inner_loop(2,2)
			do i = flow_inner_loop(1,1),flow_inner_loop(1,2)
				do dim =  1, 3
				do dim1 = 1, dimensions
				do dim2 = 1, dimensions
				!	if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim1,1),j-I_m(dim1,2),k-I_m(dim1,3)) == 0).or.(bc%bc_markers(i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)) == 0)) then 
						if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
							if(((dim1-dim) == 1).or.((dim1-dim) == -2)) then
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) - v_f%pr(dim1)%cells(dim1,i,j,k) + v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))
							else
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) + v_f%pr(dim1)%cells(dim1,i,j,k) - v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))
							end if
						end if
				!	end if
				end do
				end do

					vorticity(dim,i,j,k) = vorticity(dim,i,j,k) / cell_size(1) 
					if(abs(vorticity(dim,i,j,k)) < 1e-03) vorticity(dim,i,j,k) = 0.0_dkind
				end do
				
				
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(guided)	
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number1	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							do dim2 = 1, dimensions
								if(dim2 /= dim) then
									bound_number2	= bc%bc_markers(i+sign*I_m(dim,1)-I_m(dim2,1),j+sign*I_m(dim,2)-I_m(dim2,2),k+sign*I_m(dim,3)-I_m(dim2,3))
									bound_number3	= bc%bc_markers(i+sign*I_m(dim,1)+I_m(dim2,1),j+sign*I_m(dim,2)+I_m(dim2,2),k+sign*I_m(dim,3)+I_m(dim2,3))
									exit
								else
									bound_number2	= 0
									bound_number3	= 0
								end if
							end do
							if ((( bound_number1 /= 0 ).and.(bound_number2 /= 0)).or.(( bound_number1 /= 0).and.(bound_number3/= 0))) then
							
								boundary_type_name1 = bc%boundary_types(bound_number1)%get_type_name()
								if (bound_number2 /= 0) boundary_type_name2 = bc%boundary_types(bound_number2)%get_type_name()
								if (bound_number3 /= 0) boundary_type_name3 = bc%boundary_types(bound_number3)%get_type_name()
   
								boundary_type_name = 'outlet'
								
								if (( bound_number1 /= 0 ).and.(bound_number2 /= 0)) then
									if ((boundary_type_name1=='wall').or.(boundary_type_name2=='wall')) boundary_type_name = 'wall'
									if ((boundary_type_name1=='inlet').or.(boundary_type_name2=='inlet')) boundary_type_name = 'inlet'
								
								
									select case(boundary_type_name)
										case('wall')
											!if (sign < 0) then
											!	vorticity(3,i,j,k) = vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) !0.0_dkind
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = vorticity(3,i,j,k) !0.0_dkind
											!end if
										case('inlet')
											!if (sign < 0) then
											!	vorticity(3,i,j,k) = 0.0_dkind
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dkind
											!end if	
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i,j,k) = 0.0_dkind !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)) !0.5_dkind*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))) !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) =  cell_size(1)*(vorticity(3,i,j,k) - (v_f%pr(dim)%cells(dim,i,j,k) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) = -cell_size(1)*(vorticity(3,i,j,k) + (v_f%pr(dim)%cells(dim,i,j,k) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1)) 
												end if
											else
												vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dkind !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) !0.5_dkind*(vorticity(3,i,j,k) + vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = -cell_size(1)*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1)+I_m(dim,1),j-I_m(dim2,2)+I_m(dim,2),k-I_m(dim2,3)+I_m(dim,3)))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) =  cell_size(1)*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1)+I_m(dim,1),j-I_m(dim2,2)+I_m(dim,2),k-I_m(dim2,3)+I_m(dim,3)))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1)) 
												end if
											end if									
									end select
								end if
								
								if (( bound_number1 /= 0 ).and.(bound_number3 /= 0)) then
									if ((boundary_type_name1=='wall').or.(boundary_type_name3=='wall')) boundary_type_name = 'wall'
									if ((boundary_type_name1=='inlet').or.(boundary_type_name3=='inlet')) boundary_type_name = 'inlet'
   
									select case(boundary_type_name)
										case('wall')
											!if (sign < 0) then
											!	vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) = -vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) !0.0_dkind
											!else
											!	vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = -vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) !0.0_dkind
											!end if
										case('inlet')
											!if (sign < 0) then
											!	vorticity(3,i,j,k) = 0.0_dkind
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dkind
											!end if	
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) = 0.0_dkind !2.0*vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) - vorticity(3,i+2*I_m(dim,1)+I_m(dim2,1),j+2*I_m(dim,2)+I_m(dim2,2),k+2*I_m(dim,3)+I_m(dim2,3)) !vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))!2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)) !0.5_dkind*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))) !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) =  cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) = -cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) + (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1)) 
												end if
											else
												vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = 0.0_dkind !2.0*vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - vorticity(3,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) !vorticity(3,i,j,k)!2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) !0.5_dkind*(vorticity(3,i,j,k) + vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = -cell_size(1)*(vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) - (v_f%pr(dim)%cells(dim,i+I_m(dim2,1)+I_m(dim,1),j+I_m(dim2,2)+I_m(dim,2),k+I_m(dim2,3)+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) =  cell_size(1)*(vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) + (v_f%pr(dim)%cells(dim,i+I_m(dim2,1)+I_m(dim,1),j+I_m(dim2,2)+I_m(dim,2),k+I_m(dim2,3)+I_m(dim,3)) - v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1)) 
												end if
											end if									
									end select
								end if
							end if
						end do
					end do	
				end if
			end do
			end do
			end do
			!$omp end do

			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))
	
				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	
	
				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))
				
				!$omp do collapse(3) schedule(guided)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then 
						if (predictor) then
							if(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) > 1e-010) then
								F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* g(dim))
							end if

							if (this%viscosity_flag)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(0.5_dkind*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
						else
							if(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) > 1e-010) then
								F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) *(this%rho_0 - rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* g(dim))
							end if
								
							if (this%viscosity_flag)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) *(0.5_dkind*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
						
						end if
					end if
				end do
				end do
				end do
				!$omp end do
			end do
			
			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))
	
				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	
	
				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))
				
				!$omp do collapse(3) schedule(guided)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)
					if (j==65) then
						continue
					end if
				
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then 
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
!									F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)	![m/s^2]						!# Lagrangian droplets solver
								if (bc%bc_markers(i,j,k) == 0) then
									F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dkind*(v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum droplets solver
								end if
							end do		
						end if
							
						if (this%additional_particles_phases_number /= 0) then
							do particles_phase_counter = 1, this%additional_particles_phases_number
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)	![m/s^2]						!# Lagrangian particles solver
								!if (bc%bc_markers(i,j,k) == 0) then
								!	F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dkind*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum particles solver
								!end if
							end do		
						end if							
							
						do dim1 = 1, 3
						do dim2 = 1, dimensions
						!	if((i>10).and.(j>10).and.(i<cons_inner_loop(1,2)-10).and.(j<cons_inner_loop(2,2)-10)) then
								if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
								if(((dim-dim1) == 1).or.((dim-dim1) == -2)) then
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								else
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								end if
								end if
						!	end if
						end do
						end do
					end if
				end do
				end do
				end do
				!$omp end do
			end do
	
			!!$omp do collapse(3) schedule(guided)	
			!do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			!do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			!do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
			!	if(bc%bc_markers(i,j,k) == 0) then
			!		do dim = 1, dimensions
			!			do plus = 1,2
			!				sign			= (-1)**plus
			!				bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
			!				if( bound_number /= 0 ) then
			!					boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
			!					select case(boundary_type_name)
			!						case('outlet')
			!							!if (sign < 0) then
			!							!	F_a%cells(dim,i,j,k) = F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
			!							!else
			!							!	F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = F_a%cells(dim,i,j,k)
			!							!end if
			!							!do dim2 = 1, dimensions
			!							!	if (dim /= dim2) then
			!							!		F_a%cells(dim2,i,j,k) = F_a%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
			!							!	end if
			!							!end do
			!					end select
			!				end if
			!			end do
			!		end do	
			!	end if
			!end do
			!end do
			!end do
			!!$omp end do

			!$omp do collapse(3) schedule(guided) reduction(+:grad_F_a_summ)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					
					lame_coeffs		= 1.0_dkind				
			
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dkind
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
					end select						
					
					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)	= (F_a%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_a%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
						grad_F_a_summ		= grad_F_a_summ + grad_F_a(dim,i,j,k) 
					end do

				end if
			end do
			end do
			end do		
			!$omp end do			
			
			if (this%all_Neumann_flag) then	
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)  = grad_F_a(dim,i,j,k)  - grad_F_a_summ/this%sub_cells_number(0)
					end do
				end if
			end do
			end do
			end do		
			!$omp end do			
			end if

			!$omp end parallel
			
			overall_poisson_iteration = 0			
			pressure_iteration	= 0
			pressure_converged	= .false.
			
			do while ((.not.pressure_converged).and.(pressure_iteration < 200)) 

				this%p_old	= p_dyn%cells
				
				H_max_old	= 10.0
				a_norm_init = 0.0_dkind

				grad_F_b_summ = 0.0_dkind
				
				!$omp parallel default(none)  private(i,j,k,dim,dim2,loop,plus,sign,bound_number,boundary_type_name,residual,lame_coeffs,farfield_velocity) , &
				!$omp& firstprivate(this) , &
				!$omp& shared(F_a,F_b,grad_F_b,grad_F_b_summ,v,v_f,v_f_old,p_dyn,H,H_old,R,ddiv_v_dt,rho_old,rho_int,bc,predictor,cons_utter_loop,cons_inner_loop,dimensions,cell_size,coordinate_system,mesh,a_norm_init,time_step,farfield_velocity_array)
				
				do dim = 1, dimensions
					loop(3,1) = cons_inner_loop(3,1)
					loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

					loop(2,1) = cons_inner_loop(2,1)
					loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

					loop(1,1) = cons_inner_loop(1,1)
					loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))	

					!$omp do collapse(3) schedule(guided)
					do k = loop(3,1),loop(3,2)
					do j = loop(2,1),loop(2,2)
					do i = loop(1,1),loop(1,2)

						if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
							if (predictor) then	
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_old%cells(i,j,k))		&
														/	(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_old%cells(i,j,k))	&
														*	(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	!/ cell_size(1)
							else
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_int%cells(i,j,k))		&
														/	(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_int%cells(i,j,k))	&
														*	(1.0_dkind/rho_int%cells(i,j,k)	- 1.0_dkind/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	!/ cell_size(1)
							end if

						end if
					end do
					end do
					end do
					!$omp end do

				end do
				
				!$omp do collapse(3) schedule(guided) reduction(+:grad_F_b_summ)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dkind				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select						
					
						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)	= (F_b%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_b%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
							grad_F_b_summ	= grad_F_b_summ + grad_F_b(dim,i,j,k) 
						end do
						
					end if
				end do
				end do
				end do		
				!$omp end do			
			
				if (this%all_Neumann_flag) then	
				!$omp do collapse(3) schedule(guided)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)  = grad_F_b(dim,i,j,k)  - grad_F_b_summ/this%sub_cells_number(0)
						end do
					end if
				end do
				end do
				end do		
				!$omp end do				
				end if
				
				!$omp do collapse(3) schedule(guided) reduction(+:a_norm_init)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dkind				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select						
					
						R%cells(i,j,k)	= 0.0_dkind
						
						R%cells(i,j,k)	= R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

						do dim = 1, dimensions
								
							R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
							R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

							R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k) 

							R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k)

							do plus = 1,2
								sign			= (-1)**plus
								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
											farfield_velocity = 0.0_dkind
											if(predictor) then
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																							
											else
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)

											end if													
										case('inlet')
										!	farfield_velocity = farfield_velocity_array(factor * j)
											farfield_velocity = farfield_velocity_array(1)
											if(predictor) then
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																				- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																							
											else
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																				- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																						+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dkind * cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											end if
										case('outlet')
											R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
											if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dkind) then
												do dim2 = 1, dimensions
													R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
												end do
											end if
									end select
								end if
							end do
						end do

						a_norm_init = a_norm_init + abs(R%cells(i,j,k))
						
					end if
				end do
				end do
				end do
				!$omp end do
				!$omp end parallel	
				
				poisson_iteration	= 0 
				converged			= .false.
				H_max				= 0.0_dkind
				a_norm_prev			= a_norm_init
				beta				= 2.0_dkind/3.0_dkind
				v_cycle_converged	= .false.
			

				nu_0 = 5
				nu_1 = 2
				nu_2 = 2
				V_cycle_depth = this%number_of_meshes - 2

				
				tolerance = 1e-02_dkind
					
				v_cycle_iteration = 0
				do while (((v_cycle_iteration <= 0).or.(.not.v_cycle_converged)).and.(v_cycle_iteration <= nu_0))
			
					nu_1 = 2!(v_cycle_iteration+2) * 2
					nu_2 = 1!(v_cycle_iteration+2) * 2
				
					!$omp parallel default(none)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) ,			&
					!$omp& firstprivate(this) , &
					!$omp& shared(H_summ,H,H_old,R,a_norm,a_norm_init,a_norm_prev,H_max,H_max_old,H_average,F_a,F_b,p_dyn,rho_old,rho_int,beta,ddiv_v_dt,v,				&
					!$omp& v_f,v_f_old,cons_inner_loop,cons_utter_loop,subgrid_cons_inner_loop,bc,dimensions,cell_size,coordinate_system,mesh,converged,predictor,time_step,time,poisson_iteration,	&
					!$omp& farfield_velocity_array,tolerance,nu_1,sub_F, sub_F_summ,r_summ)
					
					!$omp master
					poisson_iteration = 0
					!$omp end master
					!$omp barrier
					
					do while ((.not.converged).and.(poisson_iteration <= nu_1))
					
						!$omp barrier
				
						!$omp master
						a_norm	= 0.0_dkind
						converged = .false.
						!$omp end master
					
						!$omp barrier

						!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
			  
								lame_coeffs		= 1.0_dkind				
		
								select case(coordinate_system)
									case ('cartesian')	
										lame_coeffs			= 1.0_dkind
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
								end select								

								R%cells(i,j,k) = 0.0_dkind
								
								R%cells(i,j,k) = R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)								

								do dim = 1, dimensions
								
									R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
									R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k) 

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k) 

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
													farfield_velocity = 0.0_dkind
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																							
													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																								
																																								
													end if													
												case('inlet')
												!	farfield_velocity = farfield_velocity_array(factor * j)
													farfield_velocity = farfield_velocity_array(1)
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																						- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																							
													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																						- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																								+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dkind * cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
													R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
													if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dkind) then
														do dim2 = 1, dimensions
															R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
														end do
													end if										
											end select		
										end if
									end do
								end do								
								
								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*R%cells(i,j,k)
  
								a_norm = a_norm + abs(R%cells(i,j,k))
								
							end if
						end do
						end do
						end do
						!$omp end do
					
						!$omp do collapse(3) schedule(guided) 
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							R%cells(i,j,k)	= R%cells(i,j,k) / cell_size(1) / cell_size(1)
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do
						
						!$omp barrier
						!$omp master
						!if(a_norm/a_norm_init < tolerance) converged = .true.
					
						if ((poisson_iteration == 0).or.(poisson_iteration == nu_1)) then
						!	print *, "T", a_norm, poisson_iteration
						end if
						!pause	
				
						poisson_iteration	= poisson_iteration + 1
						!$omp end master
						!$omp barrier
					end do

					!$omp master
					a_norm		= 0.0_dkind
					converged   = .false.
                    sub_F_summ	= 0.0_dkind
                    r_summ		= 0.0_dkind
					!$omp end master					
					
					!$omp do collapse(3) schedule(guided) reduction(+:sub_F_summ,r_summ)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						sub_F(0)%cells(i,j,k)	= R%cells(i,j,k)
                        
                        select case(coordinate_system)
							case ('cartesian')	
								sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k)
								r_summ = r_summ + 1.0_dkind! this%sub_cells_number(mesh_iter+1)
							case ('cylindrical')
								! x -> r, y -> z
								sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k) * mesh%mesh(1,i,j,k)
								r_summ		= r_summ + mesh%mesh(1,i,j,k)
							case ('spherical')
								! x -> r
								sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k) * mesh%mesh(1,i,j,k)**2
								r_summ		= r_summ + mesh%mesh(1,i,j,k)**2
						end select
                      
					end do
					end do
                    end do
					!$omp end do					

                    if (this%all_Neumann_flag) then	
					!$omp do collapse(3) schedule(guided) 
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						sub_F(0)%cells(i,j,k) = sub_F(0)%cells(i,j,k) - sub_F_summ / r_summ !this%sub_cells_number(0)
					end do
					end do
                    end do
					!$omp end do    
                    end if
                    
					!$omp end parallel
					
					call this%V_cycle(0,time_step,nu_1,nu_2,tolerance,predictor,V_cycle_depth)
		
					!$omp parallel default(none)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) ,			&
					!$omp& firstprivate(this) , &
					!$omp& shared(H_summ,H,H_old,R,a_norm,a_norm_init,a_norm_prev,H_max,H_max_old,H_average,F_a,F_b,p_dyn,rho_old,rho_int,beta,ddiv_v_dt,v,					&
					!$omp& v_f,v_f_old,cons_inner_loop,cons_utter_loop,subgrid_cons_inner_loop,bc,dimensions,cell_size,coordinate_system,mesh,converged,predictor,time_step,time,poisson_iteration,	&
					!$omp& farfield_velocity_array,tolerance,nu_2,sub_E)
					
					!$omp do collapse(3) schedule(guided)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						H_old%cells(i,j,k) = H%cells(i,j,k) + sub_E(0)%cells(i,j,k)
					end do
					end do
					end do
					!$omp end do					

					!$omp master
					a_norm	= 0.0_dkind
					!$omp end master					
					
					!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						if(bc%bc_markers(i,j,k) == 0) then
			  
							lame_coeffs		= 1.0_dkind				
		
							select case(coordinate_system)
								case ('cartesian')	
									lame_coeffs			= 1.0_dkind
								case ('cylindrical')
									! x -> r, y -> z
									lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
									lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
									lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
								case ('spherical')
									! x -> r
									lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
									lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
									lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
							end select								

							R%cells(i,j,k) = 0.0_dkind
								
							R%cells(i,j,k) = R%cells(i,j,k) +  ddiv_v_dt%cells(i,j,k) * cell_size(1) * cell_size(1)							

							do dim = 1, dimensions
								
								R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2) 
								
								R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2) 

								R%cells(i,j,k)	= R%cells(i,j,k) +	grad_F_a(dim,i,j,k) * cell_size(1) 

								R%cells(i,j,k)	= R%cells(i,j,k) +	grad_F_b(dim,i,j,k) * cell_size(1) 

								do plus = 1,2
									sign			= (-1)**plus
									bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
									if( bound_number /= 0 ) then
										boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
										select case(boundary_type_name)
											case('wall')
												farfield_velocity = 0.0_dkind
												if(predictor) then
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																							
												else
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																								
												end if													
											case('inlet')
											!	farfield_velocity = farfield_velocity_array(factor * j)
												farfield_velocity = farfield_velocity_array(1)
												if(predictor) then
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																					- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
												else
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																					- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																							+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dkind * cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
												end if
											case('outlet')
												R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
												if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dkind) then
													do dim2 = 1, dimensions
														R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
													end do
												end if
										end select		
									end if
								end do
							end do	
							
							a_norm = a_norm + abs(R%cells(i,j,k))
							
						end if
					end do
					end do
					end do			
					!$omp end do					
					
					!$omp barrier
					!$omp master
					poisson_iteration = 0
					!print *, "01", a_norm, poisson_iteration
					!$omp end master
					!$omp barrier
					
					do while ((.not.converged).and.(poisson_iteration <= nu_2))
					
						!$omp barrier
				
						!$omp master
						a_norm	= 0.0_dkind
						converged = .false.
						!$omp end master
					
						!$omp barrier

						!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
			  
								lame_coeffs		= 1.0_dkind				
		
								select case(coordinate_system)
									case ('cartesian')	
										lame_coeffs			= 1.0_dkind
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
								end select								

								R%cells(i,j,k) = 0.0_dkind
								
								R%cells(i,j,k) = R%cells(i,j,k) + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

								do dim = 1, dimensions
								
									R%cells(i,j,k)	= R%cells(i,j,k) +	(H_old%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H_old%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
									R%cells(i,j,k)	= R%cells(i,j,k) -	(H_old%cells(i,j,k) - H_old%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_a(dim,i,j,k) 

									R%cells(i,j,k)	= R%cells(i,j,k) +	cell_size(1)*grad_F_b(dim,i,j,k) 

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
													farfield_velocity = 0.0_dkind
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	

													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	

													end if													
												case('inlet')
												!	farfield_velocity = farfield_velocity_array(factor * j)
													farfield_velocity = farfield_velocity_array(1)
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																						- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)				&
																						- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																								+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dkind * cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
													R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
													if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dkind) then
														do dim2 = 1, dimensions
															R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
														end do
													end if
											end select		
										end if
									end do
								end do
								
								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*R%cells(i,j,k)
  
								a_norm = a_norm + abs(R%cells(i,j,k))
	 
								R%cells(i,j,k)	= R%cells(i,j,k) / cell_size(1) / cell_size(1)
								
							end if
						end do
						end do
						end do
						!$omp end do
					
						!$omp do collapse(3) schedule(guided) 
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do
						
						!$omp barrier
						!$omp master
						!if(a_norm/a_norm_init < tolerance) converged = .true.
					
						if ((poisson_iteration == 0).or.(poisson_iteration == nu_2)) then
						!	print *, "T", a_norm, poisson_iteration
						end if
						!pause	
				
						poisson_iteration	= poisson_iteration + 1
						!$omp end master
						!$omp barrier
					end do				
					
					!$omp end parallel						
					
					!print *, v_cycle_iteration
					v_cycle_iteration = v_cycle_iteration + 1
				
					do mesh_iter = 0, this%number_of_meshes-1
						this%sub_E_old(mesh_iter)%cells	= 0.0_dkind
						this%sub_E(mesh_iter)%cells		= 0.0_dkind
					end do	
					
					if(abs(a_norm_init) > 1e-10_dkind) then
						print *,a_norm,a_norm_init,a_norm/a_norm_init
						if((a_norm/a_norm_init < tolerance).or.(a_norm < 100e-0)) then
							v_cycle_converged = .true.
							
						end if	
					end if

					!pause
				end do
				
				!!$omp end parallel
				
				print *, 'Poisson iteration:', poisson_iteration

				overall_poisson_iteration = overall_poisson_iteration + poisson_iteration

				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
	
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1,dimensions															 
							do plus = 1,2
								sign			= (-1)**plus
								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
											if(predictor) then
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  H%cells(i,j,k) &
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)		
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  H%cells(i,j,k) &
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)
											end if
										case('outlet')
											if (sign == 1) then		!# Правая граница
												if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dkind) then		!# Выток, берутся значения слева
												
													!if (predictor) then	
													!	H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/rho_old%cells(i,j,k)
													!else
													!	H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/rho_int%cells(i,j,k)
													!end if
													!
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	- H%cells(i,j,k) 
													do dim2 = 1, dimensions
														H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + v%pr(dim2)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))**2
													end do
												else																					!# Вток, берутся значения справа
													farfield_density		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_density()
													farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													farfield_pressure		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
							
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/farfield_density + 0.5_dkind*(farfield_velocity **2)
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	- H%cells(i,j,k) 
													!do dim2 = 1, dimensions
													!	H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + v%pr(dim2)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))**2
													!end do
													
												end if
											else					!# Левая граница
												if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dkind) then										!# Выток, берутся значения справа
										 
													!if (predictor) then	
													!	H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/rho_old%cells(i,j,k)
													!else
													!	H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/rho_int%cells(i,j,k)
													!end if
													!												
													H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	- H%cells(i,j,k) 
													do dim2 = 1, dimensions
														H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) = H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v%pr(dim2)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2
													end do
												else																					!# Вток, берутся значения справа
													farfield_density		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_density()
													farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													farfield_pressure		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
												
													H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/farfield_density + 0.5_dkind*(farfield_velocity **2) 
													H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	- H%cells(i,j,k) 
													!do dim2 = 1, dimensions
													!	H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) = H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v%pr(dim2)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2
													!end do
													
												end if
											end if
												
										case('inlet')
	 
										!	farfield_velocity = farfield_velocity_array(factor * j)
											farfield_velocity = farfield_velocity_array(1)
												
											if (predictor) then
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i,j,k)														&
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)	&	
																												- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)/time_step
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i,j,k)														&
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dkind * cell_size(1)	&	
																												- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																														+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dkind * cell_size(1)/(0.5_dkind*time_step)
											end if
									end select
								end if
							end do	
						end do
					end if
				end do
				end do
				end do	
				
				call this%calculate_dynamic_pressure(time_step,predictor)
				pressure_converged = .true.				
				
				!$omp parallel default(none)  private(i,j,k,dim,H_residual,r_i,r_j,r_k,plus,sign,bound_number,boundary_type_name,farfield_pressure) , &
				!$omp& firstprivate(this) , &
				!$omp& shared(H_summ,p_dyn,rho_old,rho_int,H,v_f,v_f_old,F_a,F_b,predictor,pressure_converged,cons_inner_loop,cons_utter_loop,bc,dimensions,cell_size,time_step,time)

				
				!$omp do collapse(3) schedule(guided)	reduction(.and.:pressure_converged)	
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1,dimensions
							if (predictor) then
								if (pressure_converged) then
								if (abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>20.0_dkind/cell_size(1)/cell_size(1)) then
									pressure_converged = .false.
									print *, 'Pressure error', abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))
								end if
								end if
							else
								if (pressure_converged) then
								if (abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_int%cells(i,j,k)	- 1.0_dkind/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>20.0_dkind/cell_size(1)/cell_size(1)) then
									pressure_converged = .false.
									print *, 'Pressure error',  abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_int%cells(i,j,k)	- 1.0_dkind/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))
								end if
								end if						
							end if
						end do
	
					end if
				end do
				end do
				end do		
				!$omp end do
				!$omp end parallel

				pressure_iteration = pressure_iteration + 1

				continue
			end do
				
			print *, 'Pressure iteration:', pressure_iteration
			print *, 'Overall pressure iteration:', overall_poisson_iteration
			
			continue
			
		end associate

	end subroutine
	
	subroutine calculate_dynamic_pressure(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor	
				
		real(dkind)	:: farfield_density, farfield_pressure, farfield_velocity
		real(dkind)	:: vel_abs
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions 
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()

		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
	
		associate (	rho_old			=> this%rho_old%s_ptr		, &
					rho_int			=> this%rho_int%s_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					H				=> this%H%s_ptr				, &
					v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					bc				=> this%boundary%bc_ptr)
	
			!$omp parallel default(none)  private(i,j,k,dim,vel_abs,plus,sign,bound_number,boundary_type_name,farfield_density,farfield_pressure) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(p_dyn,v,v_f,rho_old,rho_int,H,predictor,cons_inner_loop,bc,dimensions,cell_size)	
					
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then

					vel_abs				= 0.0_dkind
					do dim = 1,dimensions
						vel_abs				= vel_abs + (0.5_dkind*(v_f%pr(dim)%cells(dim,i,j,k) + v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))))**2	
					end do	
					
					if (predictor) then	
						p_dyn%cells(i,j,k)	=  (H%cells(i,j,k) - 0.5_dkind * vel_abs)*rho_old%cells(i,j,k)
					else
						p_dyn%cells(i,j,k)	=  (H%cells(i,j,k) - 0.5_dkind * vel_abs)*rho_int%cells(i,j,k)
					end if
	
					do dim = 1,dimensions															 
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dkind*(	v%pr(dim)%cells(i,j,k) **2)  
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if
                                    case('outlet')
										!if (sign == 1) then		!# Правая граница
										!	if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dkind) then		!# Выток, берутся значения слева
										!		p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	else																					!# Вток, берутся значения справа
										!		p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	end if
										!else					!# Левая граница
										!	if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dkind) then										!# Выток, берутся значения справа
										!		p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	else																					!# Вток, берутся значения справа
										!		p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	end if
										!end if
											
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dkind*(	v%pr(dim)%cells(i,j,k) **2)   
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if                                        
                                        
									case('inlet')
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dkind*(	v%pr(dim)%cells(i,j,k) **2)   
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if
								end select
							end if
						end do	
					end do
				end if
			end do
			end do
			end do	
			
			!$omp end do
			!$omp end parallel					
			
		end associate
	end subroutine
	
	subroutine calculate_velocity(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor

		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions, iterations
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		integer					:: plus, sign, bound_number
		character(len=20)		:: boundary_type_name
		
		integer :: i,j,k,dim,dim1,dim2,spec
	
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate (	v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					H				=> this%H%s_ptr				, &					
					bc				=> this%boundary%bc_ptr)

			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))		
				
				do k = loop(3,1), loop(3,2)		
				do j = loop(2,1), loop(2,2)		
				do i = loop(1,1), loop(1,2)		
				
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

						if (predictor)	then
							v_f_old%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k)
							v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) - time_step * (F_a%cells(dim,i,j,k) + F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(1))
						else
							v_f%pr(dim)%cells(dim,i,j,k) = 0.5_dkind * (v_f%pr(dim)%cells(dim,i,j,k) + v_f_old%pr(dim)%cells(dim,i,j,k)) - (0.5_dkind * time_step) * (F_a%cells(dim,i,j,k) + F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(1) )
						end if	

					end if

				end do
				end do
				end do
			end do
				
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')
									!	v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) = 0.0_dkind 
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))
									case('outlet')
										!if (sign > 0 ) then
										!	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f%pr(dim)%cells(dim,i,j,k)
										!else
										!	v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) 
										!end if
										
										!v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))
										
										!do dim2 = 1, dimensions
										!	if (dim2 /= dim) then
										!		v_f%pr(dim2)%cells(dim2,i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = 2.0_dkind*v_f%pr(dim2)%cells(dim2,i,j,k) - v_f%pr(dim2)%cells(dim2,i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))
										!	end if
										!end do
								end select
							end if
						end do
						v%pr(dim)%cells(i,j,k) = 0.5_dkind * ( v_f%pr(dim)%cells(dim,i,j,k) + v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
					end do	
				end if
				
			end do
			end do
			end do
			
		!	print *, '4', v_f%pr(2)%cells(2,50,20,1), v_f%pr(1)%cells(1,50,20,1)
			
			continue
		end associate
	end subroutine
	
	subroutine perturb_velocity_field(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step

		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions, iterations
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		real(dkind)	:: kappa_turb, gamma_turb, kx_turb, ky_turb	
		
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate (	v_f				=> this%v_f%v_ptr			, &
					bc				=> this%boundary%bc_ptr)

			kappa_turb	= 350.0e00_dkind

			call RANDOM_SEED()
			call RANDOM_NUMBER(gamma_turb)
			
		!	gamma_turb = 0.25
			gamma_turb = 2.0_dkind * pi * gamma_turb   !2.0_dkind*gamma_turb - 1.0_dkind   !cos(alpha) sin(alpha)=sqrt(1-gamma_turb**2.0)
			
			
			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

				loop(2,1) = cons_inner_loop(2,1)+5
				loop(2,2) = (cons_utter_loop(2,2)-5)*I_m(dim,2) + (cons_inner_loop(2,2)-5)*(1 - I_m(dim,2))	

				loop(1,1) = cons_inner_loop(1,1)+5
				loop(1,2) = (cons_utter_loop(1,2)-5)*I_m(dim,1) + (cons_inner_loop(1,2)-5)*(1 - I_m(dim,1))		
				
				do k = loop(3,1), loop(3,2)		
				do j = loop(2,1), loop(2,2)		
				do i = loop(1,1), loop(1,2)		
				!	if((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then	

						!call RANDOM_NUMBER(gamma_turb)
						!gamma_turb = 2.0_dkind*gamma_turb - 1.0_dkind
						!v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) + kappa_turb*gamma_turb*sqrt(this%time_step)
							
						kx_turb = 2.0_dkind*pi*sin(gamma_turb)/0.005_dkind    !0.002   !sqrt(1.0_dkind-gamma_turb**2.0)/0.002  
						ky_turb = 2.0_dkind*pi*cos(gamma_turb)/0.005_dkind   !0.002   !gamma_turb/0.002
							
						if(dim == 1)then
							v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) + kappa_turb*cos(gamma_turb)*sqrt(this%time_step)*cos(kx_turb*(i-1.0)*cell_size(1)+ky_turb*(j-0.5)*cell_size(1))
						endif
						if(dim == 2)then
							v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) - kappa_turb*sin(gamma_turb)*sqrt(this%time_step)*cos(kx_turb*(i-0.5)*cell_size(1)+ky_turb*(j-1.0)*cell_size(1))
						endif                        
				!	end if
				end do
				end do
				end do
			end do

			continue
		end associate
	end subroutine
	
	
	subroutine calculate_interm_Y_corrector(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		
		real(dkind)	:: rhs, B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dkind)	:: spec_summ
		
		real(dkind), dimension (3,3):: lame_coeffs
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)	:: coordinate_system
		
		integer	:: bound_number
		integer	:: i,j,k,dim,spec,droplets_phase_counter

		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		associate (	rho				=> this%rho%s_ptr			, &
					rho_int			=> this%rho_int%s_ptr		, &
					v_f				=> this%v_f%v_ptr			, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_droplets	=> this%Y_prod_droplets	, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

					
			!$omp parallel default(none)  private(i,j,k,dim,spec,spec_summ,rhs,flux_right,flux_left,lame_coeffs) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(v_f,Y_prod_diff,Y_prod_chem,Y,Y_old,Y_int,rho,rho_old,rho_int,bc,cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,mesh,time_step)					
			!$omp do collapse(3) schedule(guided)			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
					do spec = 1, species_number

						lame_coeffs		= 1.0_dkind				
				  
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select						

						Y_rho_int(spec) = 0.5_dkind * ( rho_old%cells(i,j,k)*Y_old%pr(spec)%cells(i,j,k) + rho_int%cells(i,j,k)*Y_int%pr(spec)%cells(i,j,k)) 
						
						rhs = 0.0_dkind
						do dim = 1, dimensions	
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho_int%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y_int%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
									flux_right = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
								else
									flux_right = rho_int%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * Y_int%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end if
							end if
						
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
								flux_left = this%CHARM_flux_limiter(rho_int%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3))*		&
																	Y_int%pr(spec)%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3)),	&
																	v_f%pr(dim)%cells(dim,i,j,k))
							else
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dkind) then
									flux_left = rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * Y_int%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								else
									flux_left = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
								end if
							end if					

							rhs = rhs - (	flux_right * v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3)&
										-	flux_left * v_f%pr(dim)%cells(dim,i,j,k)* lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2)
							
							continue
						end do	
						
						Y_rho_int(spec) = Y_rho_int(spec) + rhs * (0.5_dkind * time_step) 
						
						if (this%diffusion_flag)	Y_rho_int(spec) = Y_rho_int(spec) + 0.5_dkind * Y_prod_diff%pr(spec)%cells(i,j,k) * time_step
						
						if (this%reactive_flag)		Y_rho_int(spec) = Y_rho_int(spec) + 0.5_dkind * Y_prod_chem%pr(spec)%cells(i,j,k) * time_step
						
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
								Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_droplets(droplets_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step
							end do		
						end if									
						
					end do	
					
					rho%cells(i,j,k) = sum(Y_rho_int)

					do spec = 1, species_number
						Y%pr(spec)%cells(i,j,k) = Y_rho_int(spec) / rho%cells(i,j,k)
					end do
					
				end if
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(guided)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					spec_summ = 0.0_dkind
					do spec = 1,species_number
						spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dkind)
					end do
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ
					end do
				end if
			end do
			end do
			end do
			!$omp end do
			!$omp end parallel
			
		!	print *, '5', Y%pr(1)%cells(50,20,1),rho%cells(50,20,1)
			
		end associate

    end subroutine	

	
	subroutine if_stabilized(this,time,stabilized)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time    

		logical				,intent(out)	:: stabilized
		
        logical						:: boundary 
		real(dkind)	,dimension(3)	:: cell_size		
		
		real(dkind)					:: max_val, left_val, right_val, flame_velocity, flame_surface_length, surface_factor
		real(dkind)					:: a, b 
		real(dkind)					:: time_diff, time_delay, time_stabilization
		real(dkind), save			:: previous_flame_location = 0.0_dkind, current_flame_location = 0.0_dkind, farfield_velocity = 0.0_dkind
		real(dkind), save			:: previous_time = 0.0_dkind, current_time = 0.0_dkind
        real(dkind), save			:: av_flame_velocity = 0.0_dkind, previous_av_flame_velocity = 0.0_dkind
		real(dkind), dimension(20),	save	:: flame_velocity_array = 0.0_dkind
		integer		,save			:: correction = 0, counter = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dkind)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2 
		
		integer :: CO_index, H2O2_index, HO2_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

!		CO_index	= this%chem%chem_ptr%get_chemical_specie_index('CO')
!		HO2_index	= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		
		stabilized = .false.
		
		associate (	v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr)

	!	time_delay			= 1e-04_dkind			
	!	time_diff			= 5e-05_dkind
	!	time_stabilization	= 5e-04_dkind			
					
		time_delay			= 1e-05_dkind			
		time_diff			= 1e-04_dkind
		time_stabilization	= 5e-06_dkind			
		
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# 1D front tracer
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				max_val = 0.0_dkind
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
					
						!! Grad temp
						!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val) then
						!	max_val = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
					
						! max H
						if (abs(Y%pr(HO2_index)%cells(i,j,k)) > max_val) then
							max_val = Y%pr(HO2_index)%cells(i,j,k)
							flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
							flame_front_index = i
						end if					
					
						! max CO
						!if (abs(Y%pr(CO_index)%cells(i,j,k)) > max_val) then
						!	max_val = Y%pr(CO_index)%cells(i,j,k)
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
				
					end if
				end do
			end do
			end do
	
			!left_val	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_val	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
			
			!left_val	= Y%pr(CO_index)%cells(flame_front_index-1,1,1)
			!right_val	= Y%pr(CO_index)%cells(flame_front_index+1,1,1)	
			 
			left_val	= Y%pr(HO2_index)%cells(flame_front_index-1,1,1)
			right_val	= Y%pr(HO2_index)%cells(flame_front_index+1,1,1)				 
			
			a = (right_val + left_val - 2.0_dkind * max_val)/2.0_dkind/cell_size(1)**2
			b = (max_val - left_val)/cell_size(1) - a*(2.0_dkind*flame_front_coords(1) - cell_size(1))
			
			current_flame_location = -b/2.0_dkind/a

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

			if (boundary) then
				stabilized = .true.
			end if
			
			correction = correction + 1
			
		end if	
			
		end associate
    end subroutine	
    
    
	subroutine stabilizing_inlet_1D(this,time,stabilized)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time
		logical				,intent(out)	:: stabilized
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		real(dkind)					:: max_val, left_val, right_val, flame_velocity, flame_surface_length, surface_factor
		real(dkind)					:: a, b 
		real(dkind)					:: time_diff, time_delay, time_stabilization
		real(dkind), save			:: previous_flame_location = 0.0_dkind, current_flame_location = 0.0_dkind, farfield_velocity = 0.0_dkind
		real(dkind), save			:: previous_time = 0.0_dkind, current_time = 0.0_dkind
		integer		,save			:: correction = 0, counter = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dkind)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2 
		
		integer :: CO_index, H2O2_index, HO2_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

!		CO_index	= this%chem%chem_ptr%get_chemical_specie_index('CO')
!		HO2_index	= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		
		stabilized = .false.
		
		associate (	v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr)

	!	time_delay			= 1e-04_dkind			
	!	time_diff			= 5e-05_dkind
	!	time_stabilization	= 5e-04_dkind			
					
		time_delay			= 1e-05_dkind			
		time_diff			= 1e-05_dkind
		time_stabilization	= 5e-06_dkind			
		
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# 1D front tracer
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				max_val = 0.0_dkind
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
					
						!! Grad temp
						!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val) then
						!	max_val = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
					
						! max H
						if (abs(Y%pr(HO2_index)%cells(i,j,k)) > max_val) then
							max_val = Y%pr(HO2_index)%cells(i,j,k)
							flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
							flame_front_index = i
						end if					
					
						! max CO
						!if (abs(Y%pr(CO_index)%cells(i,j,k)) > max_val) then
						!	max_val = Y%pr(CO_index)%cells(i,j,k)
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
				
					end if
				end do
			end do
			end do
	
			!left_val	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_val	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
			
			!left_val	= Y%pr(CO_index)%cells(flame_front_index-1,1,1)
			!right_val	= Y%pr(CO_index)%cells(flame_front_index+1,1,1)	
			 
			left_val	= Y%pr(HO2_index)%cells(flame_front_index-1,1,1)
			right_val	= Y%pr(HO2_index)%cells(flame_front_index+1,1,1)				 
			
			a = (right_val + left_val - 2.0_dkind * max_val)/2.0_dkind/cell_size(1)**2
			b = (max_val - left_val)/cell_size(1) - a*(2.0_dkind*flame_front_coords(1) - cell_size(1))
			
			current_flame_location = -b/2.0_dkind/a

			flame_velocity = 0.0_dkind
			if(correction == 0) then
				farfield_velocity = farfield_velocity_array(1)
				previous_flame_location = current_flame_location
			end if
			
			if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
				flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)

		!		farfield_velocity_array = max(farfield_velocity_array - 0.25_dkind*flame_velocity,0.0_dkind)
				
				farfield_velocity_array = farfield_velocity_array - 0.75_dkind*flame_velocity
				
				
				write (flame_loc_unit,'(4E14.6)') time, current_flame_location, flame_velocity, farfield_velocity_array(1)
				
				previous_flame_location = current_flame_location
				previous_time = current_time
				
				if( (correction /= 0).and.(abs(flame_velocity) < 1e-05))then 
					counter = counter + 1
				end if
				
			end if

			if(counter > 10) then
				stabilized = .true.
			end if
			
			correction = correction + 1
			
		end if	
			
		end associate

    end subroutine	
   
	subroutine stabilizing_inlet(this,time)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		real(dkind)					:: max_grad_temp, left_grad_temp, right_grad_temp, max_CO, left_CO, right_CO, flame_velocity, flame_surface_length, surface_factor
		real(dkind)					:: a, b 
		real(dkind)					:: time_diff, time_delay
		real(dkind), save			:: previous_flame_location = 0.0_dkind, current_flame_location = 0.0_dkind, farfield_velocity = 0.0_dkind
		real(dkind), save			:: previous_time = 0.0_dkind, current_time = 0.0_dkind
		integer		,save			:: correction = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dkind)	,dimension(2000,2)	:: lp_coord
		real(dkind)	,dimension(2)		:: lp_copy	
		integer(dkind)	,dimension(7,2000)	:: chains
		integer			,dimension(7)		:: chain_length		
		
		real(dkind)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2, chain_index, lp_chain
		
		integer :: CO_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		CO_index	= this%chem%chem_ptr%get_chemical_specie_index('CO')
		
		associate (	v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					bc				=> this%boundary%bc_ptr)

		time_delay	= 1e-04_dkind			
		time_diff	= 5e-05_dkind
					
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# Simple front tracer
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				max_grad_temp = 0.0_dkind
				max_CO = 0.0_dkind
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
					
						!! Grad temp
						!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_grad_temp) then
						!	max_grad_temp = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
						!	flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
						!	flame_front_index = i
						!end if
					
						! max CO
						if (abs(Y%pr(CO_index)%cells(i,j,k)) > max_CO) then
							max_CO = Y%pr(CO_index)%cells(i,j,k)
							flame_front_coords(j) = (i - 0.5_dkind)*cell_size(1) 
							flame_front_index = i
						end if
				
					end if
				end do
			end do
			end do
			
			!# 2D front tracer
	!		T_flame = 1000.0_dkind
	!	
	!		lp_index = 0 
	!		do j = cons_inner_loop(2,1), cons_inner_loop(2,2)-1
	!		do i = cons_inner_loop(1,1)+1, cons_inner_loop(1,2)-1
 !
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	- 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	- 0.5*cell_size(1)	
	!			end if
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	- 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)	+ 0.5*cell_size(1)	
	!			end if
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	+ 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	+ 0.5*cell_size(1)		
	!			end if
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	+ 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	- 0.5*cell_size(1)		
	!			end if
 !
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	- 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	- 0.5*cell_size(1)	
	!			end if
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	- 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)	+ 0.5*cell_size(1)	
	!			end if
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	+ 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	+ 0.5*cell_size(1)		
	!			end if
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	+ 0.5*cell_size(1)
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1) 
	!				lp_coord(lp_index,2) = j*cell_size(1)	- 0.5*cell_size(1)		
	!			end if			
 !
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!			end if
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)
	!			end if
	!			if((T%cells(i,j,1) < T_flame ).and.( T%cells(i+1,j,1) < T_flame ).and.( T%cells(i,j+1,1) > T_flame ).and.(T%cells(i+1,j+1,1) > T_flame )) then				
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!			end if
	!			if((T%cells(i,j,1) > T_flame ).and.( T%cells(i+1,j,1) > T_flame ).and.( T%cells(i,j+1,1) < T_flame ).and.(T%cells(i+1,j+1,1) < T_flame )) then	
	!				lp_index = lp_index + 1
	!				lp_coord(lp_index,1) = i*cell_size(1)	
	!				lp_coord(lp_index,2) = j*cell_size(1)	
	!			end if	
	!		
	!		end do
	!		end do
	!		
	!		lp_number = lp_index
	!	
	!		lp_copies = 0
	!		do lp_index	= 1, lp_number-1 - lp_copies
	!			do lp_index2	= lp_index+1, lp_number - lp_copies
	!				if( ( abs(lp_coord(lp_index,1) - lp_coord(lp_index2,1)) < 1.0e-08_dkind)) then
	!				if( ( abs(lp_coord(lp_index,2) - lp_coord(lp_index2,2)) < 1.0e-08_dkind)) then
	!					lp_copies = lp_copies + 1
	!					do lp_index3	= lp_index2, lp_number - lp_copies
	!						lp_coord(lp_index3,:) = lp_coord(lp_index3+1,:)
	!					end do
	!				end if
	!				end if
	!			end do
	!		end do
	!	
	!		lp_number = lp_number - lp_copies
	!	
	!		lp_coord(lp_number+1:,:) = -1.0_dkind
	!		
	!		min_y = 100
	!		max_x = 0
	!		do lp_index		= 1, lp_number-1
	!			if ((lp_coord(lp_index,2) - 0.1*(lp_coord(lp_index,1) - minval(lp_coord(:lp_number,1))) < min_y)) then !.and.(lp_coord(lp_index,1) > max_x)) 
	!				min_y = lp_coord(lp_index,2) - 0.1*(lp_coord(lp_index,1) - minval(lp_coord(:lp_number,1)))
	!				lp_start = lp_index
	!			end if 
	!		end do
	!	
	!		lp_copy					= lp_coord(1,:)
	!		lp_coord(1,:)			= lp_coord(lp_start,:)
	!		lp_coord(lp_start,:)	= lp_copy		
	!	
	!!		print *, lp_number, '1' 
	!		do lp_index		= 1, lp_number-1
	!			min_dist = 100
	!			do lp_index2	= lp_index+1, lp_number
	!				lp_dist = sqrt((lp_coord(lp_index,1)-lp_coord(lp_index2,1))**2 + (lp_coord(lp_index,2)-lp_coord(lp_index2,2))**2)
	!				if (lp_dist < min_dist) then
	!					min_dist = lp_dist
	!					lp_neighbour = lp_index2
	!				end if 
	!			end do
	!			lp_copy						= lp_coord(lp_index+1,:)
	!			lp_coord(lp_index+1,:)		= lp_coord(lp_neighbour,:)
	!			lp_coord(lp_neighbour,:)	= lp_copy
	!		end do
	!	
	!		chain_index = 1
	!		lp_chain = 1
	!		chains(1,lp_chain) = 1
	!		chain_length = 1 
	!		do lp_index		= 2, lp_number-1
	!			lp_dist = sqrt((lp_coord(lp_index,1)-lp_coord(lp_index+1,1))**2 + (lp_coord(lp_index,2)-lp_coord(lp_index+1,2))**2)
	!			if (lp_dist > 2*cell_size(1)) then
	!				chain_index = chain_index + 1
	!				lp_chain = 1 
	!				chains(chain_index,lp_chain) = lp_index
	!			else
	!				lp_chain = lp_chain + 1
	!				chains(chain_index,lp_chain) = lp_index
	!				chain_length(chain_index) = chain_length(chain_index) + 1
	!			end if
	!		end do		
	!	
	!		lp_chain = chains(maxloc(chain_length,1),1)
	!		lp_coord(:lp_chain,:)						= -1.0_dkind
	!	
	!		lp_chain = chains(maxloc(chain_length,1),maxval(chain_length))
	!		lp_coord(lp_chain:,:)						= -1.0_dkind
	!		
	!		
	!!		print *, lp_number, '2'
	!		flame_surface_length = 0.0_dkind
	!		current_flame_location = 0.0_dkind
	!		lp_number2 = 0
	!		do lp_index = 1, lp_number-1
	!			if((lp_coord(lp_index,1) > 0.0_dkind).and.(lp_coord(lp_index+1,1) > 0.0_dkind)) then
	!				flame_surface_length = flame_surface_length + sqrt(	(lp_coord(lp_index,1) - lp_coord(lp_index+1,1))**2  + &
	!																	(lp_coord(lp_index,2) - lp_coord(lp_index+1,2))**2)	
	!				current_flame_location = current_flame_location + lp_coord(lp_index,1) 
	!				lp_number2 = lp_number2 + 1
	!			end if
	!		end do
	!		
	!		current_flame_location = current_flame_location/lp_number2
			
			!flame_surface_length = 0.0_dkind
			!do j = cons_inner_loop(2,1),cons_inner_loop(2,2)-1 
			!	flame_surface_length = flame_surface_length + sqrt((flame_front_coords(j)-flame_front_coords(j+1))**2 + cell_size(1)**2)
			!end do			
			
			!left_grad_temp	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_grad_temp	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
   !
			!a = (right_grad_temp + left_grad_temp - 2.0_dkind * max_grad_temp)/2.0_dkind/cell_size(1)**2
			!b = (max_grad_temp-left_grad_temp)/cell_size(1) - a*(2.0_dkind*flame_front_coords(1) - cell_size(1))			
			
			 left_CO	= Y%pr(CO_index)%cells(flame_front_index-1,1,1)
			 right_CO	= Y%pr(CO_index)%cells(flame_front_index+1,1,1)		
			
			 a = (right_CO + left_CO - 2.0_dkind * max_CO)/2.0_dkind/cell_size(1)**2
			 b = (max_CO - left_CO)/cell_size(1) - a*(2.0_dkind*flame_front_coords(1) - cell_size(1))
			
			 current_flame_location = -b/2.0_dkind/a
			
			!current_flame_location = sum(flame_front_coords) / (cons_inner_loop(2,2) - cons_inner_loop(2,1) + 1)
			!print *, lp_number, '3'
			!current_flame_location = sum(lp_coord(:lp_number,1)) / (lp_number)

			flame_velocity = 0.0_dkind
			if(correction == 0) then
				farfield_velocity = farfield_velocity_array(1)
				previous_flame_location = current_flame_location
			end if

		!	print *, correction, time
		!	pause 
			
			if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
				flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)

		!	print *, flame_velocity, farfield_velocity_array(1), current_flame_location, previous_flame_location,current_time,previous_time
		!	print *, '----'
		!	print *, flame_front_coords
		!	pause
			
				farfield_velocity_array = farfield_velocity_array - 0.25_dkind*flame_velocity
				
		!		surface_factor =  flame_surface_length / (cell_size(1) * (cons_inner_loop(2,2)- cons_inner_loop(2,1)))
		!		farfield_velocity_array = farfield_velocity * surface_factor				
				
				!write(file_name,'(A,E14.7,A)') 'front_data/front_data_',time,'.dat' 
				!open(newunit = flame_structure_unit, file = file_name, status = 'replace', form = 'formatted')
				
				!do lp_index		= 1, lp_number
				!	if(lp_coord(lp_index,1) > 0.0) then
				!		write(flame_structure_unit,'(E14.7,A,E14.7)') lp_coord(lp_index,1), '  ', lp_coord(lp_index,2)
				!	end if
				!end do
				!
				!close(flame_structure_unit)
				
				!
				!pause

				write (flame_loc_unit,'(4E14.6)') time, current_flame_location, flame_velocity, farfield_velocity_array(1)! , flame_surface_length, surface_factor
				
				previous_flame_location = current_flame_location
				previous_time = current_time
				
			end if

			correction = correction + 1
			
		end if	
			
		end associate

	end subroutine
	
	subroutine apply_boundary_conditions(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		real(dkind)					:: wall_temperature, farfield_temperature, farfield_pressure, farfield_velocity
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec
		
		character(len=20)		:: boundary_type_name
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate (	rho				=> this%rho%s_ptr			, &
					rho_int			=> this%rho_int%s_ptr		, &
					rho_old			=> this%rho_old%s_ptr		, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_old			=> this%Y_old%v_ptr			, &
					v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					h_s				=> this%h_s%s_ptr			, &
					T				=> this%T%s_ptr				, &
					bc				=> this%boundary%bc_ptr)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
					
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then
								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								
								select case(boundary_type_name)
									case('wall')
										if(this%boundary%bc_ptr%boundary_types(bound_number)%is_conductive()) then
											wall_temperature = this%boundary%bc_ptr%boundary_types(bound_number)%get_wall_temperature()
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = wall_temperature
										else
											T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = T%cells(i,j,k)
										end if	
								
										if (predictor) then
											rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_int%cells(i,j,k) * T%cells(i,j,k) / T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_old%cells(i,j,k) * T%cells(i,j,k) / T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											do spec = 1, species_number
												Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y_int%pr(spec)%cells(i,j,k)
												Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
												Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											end do
										else
											rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))					= rho%cells(i,j,k) * T%cells(i,j,k) / T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											do spec = 1, species_number
												Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= Y%pr(spec)%cells(i,j,k)
											end do
										end if
								
										h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= h_s%cells(i,j,k) * T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) / T%cells(i,j,k)

										do dim1 = 1, dimensions
											if(dim1 == dim) then
											v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= -v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											end if
										end do
										
									case('outlet')
										farfield_temperature											= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T%cells(i,j,k) !farfield_temperature
										farfield_pressure												= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
										
										if (predictor) then
											rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_int%cells(i,j,k) !rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho_old%cells(i,j,k) !rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											do spec = 1, species_number
												Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y_int%pr(spec)%cells(i,j,k) !Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
												Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y_old%pr(spec)%cells(i,j,k) !Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											end do

										else
											rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = rho%cells(i,j,k)
											do spec = 1, species_number
												Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i,j,k)
											end do
										end if
										
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											end if
										end do
										
									case('inlet')
										if(sign == 1) then
											if(v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T%cells(i,j,k)
											else
												farfield_temperature											= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= farfield_temperature
											end if
										else	
											if(v_f%pr(dim)%cells(dim,i,j,k) < 0.0_dkind) then
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T%cells(i,j,k)
											else
												farfield_temperature											= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= farfield_temperature
											end if	
										end if
										
										farfield_velocity		=  farfield_velocity_array(j)

										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  farfield_velocity
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  0.0_dkind 
											end if
										end do

										if (predictor) then
											rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											do spec = 1, species_number
												Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
												Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											end do
										end if
								end select
							end if
						end do
					end do	
				end if
			end do
			end do
			end do

			continue
			
		end associate

	end subroutine
	
	subroutine apply_poisson_boundary_conditions(this,time_step,predictor)
		class(fds_solver)	,intent(in)	:: this
		real(dkind)			,intent(in)	:: time_step
		logical				,intent(in)	:: predictor
		
		real(dkind)	:: H_residual
		real(dkind)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions 
		integer	,dimension(3,2)	:: cons_inner_loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
	
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
	
		associate (	rho_old			=> this%rho_old%s_ptr		, &
					rho_int			=> this%rho_int%s_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					H				=> this%H%s_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					bc				=> this%boundary%bc_ptr)
	
		end associate
	end subroutine

	pure function CHARM_flux_limiter(this,scalar_array,velocity)
		class(fds_solver)					,intent(in)	:: this
		real(dkind)			,dimension(4)	,intent(in)	:: scalar_array
		real(dkind)							,intent(in)	:: velocity
		real(dkind)							:: CHARM_flux_limiter
		
		real(dkind)	:: B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		
		phi_up	= 0.0_dkind
		B_r		= 0.0_dkind

		phi_loc = scalar_array(3)	- scalar_array(2)
		if (velocity > 0.0_dkind) then
			phi_up = scalar_array(2)	- scalar_array(1)
		else
			phi_up = scalar_array(4)	- scalar_array(3)
		end if
						
		if (phi_loc /= 0) then						
			r	= phi_up / phi_loc
		else
			r	= 0.0_dkind
		end if	
  
		!# CHARM
		if ( r > 1.0e-010_dkind) then
			s	= 1.0_dkind/r
			B_r	= s * (3.0_dkind*s + 1.0_dkind) / (s+1.0_dkind) / (s+1.0_dkind)	
		else
			B_r = 0.0_dkind
		end if

		!# SUPERBEE
		!B_r	= max(max(0.0_dkind,min(2.0_dkind*r,1.0_dkind)),min(r,2.0_dkind))
		
		!# Godunov
		!B_r	= 0.0_dkind
        
 		!# Central Difference
		!B_r	= 1.0_dkind       

		
		if (velocity > 0.0_dkind) then
			CHARM_flux_limiter = scalar_array(2) + B_r * 0.5_dkind * phi_up		!# CHARM
		!	CHARM_flux_limiter = scalar_array(2) + B_r * 0.5_dkind * phi_loc	!# Others
		else
			CHARM_flux_limiter = scalar_array(3) - B_r * 0.5_dkind * phi_up		!# CHARM
		!	CHARM_flux_limiter = scalar_array(3) - B_r * 0.5_dkind * phi_loc	!# Others
		end if
		
	end function	

	subroutine calculate_time_step(this)
		class(fds_solver)	,intent(inout)	:: this
		
		real(dkind)	:: delta_t_interm, delta_t_interm1, delta_t_interm1x, delta_t_interm2, delta_t_interm3, time_step, time_step2, velocity_value, divergence_value

		real(dkind)	,dimension(3)	:: cell_size

		real(dkind)	,save	:: time

		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: sign
		integer :: i,j,k,dim, spec

		time_step		= 10.0_dkind !this%initial_time_step
		time_step2		= 1000.0_dkind !	
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()			


		associate(  v				=> this%v%v_ptr				, &
					div_v_int		=> this%div_v_int%s_ptr		, &
					v_s				=> this%v_s%s_ptr			, &
					nu				=> this%nu%s_ptr			, &
					kappa			=> this%kappa%s_ptr			, &
					rho				=> this%rho%s_ptr			, &
					D				=> this%D%v_ptr				, &
					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)
		
		!$omp parallel default(none)  private(i,j,k,dim,delta_t_interm,velocity_value) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v,rho,nu,D,div_v_int,bc,time_step,cons_inner_loop,dimensions,cell_size)
		!$omp do collapse(3) schedule(guided) reduction(min:time_step)						
						
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dkind
				do dim = 1,dimensions
					velocity_value	= velocity_value + abs(v%pr(dim)%cells(i,j,k))/(cell_size(1))			!#L1 norm
				!	velocity_value	= velocity_value + (v%pr(dim)%cells(i,j,k)/minval(cell_size(1)))**2		!#L2 norm
				end do

				if((velocity_value > 0.0_dkind).or.(abs(div_v_int%cells(i,j,k)) > 0.0_dkind)) then
					delta_t_interm = 1.0_dkind/(velocity_value + abs(div_v_int%cells(i,j,k)))		!# L1 norm
				!	delta_t_interm = 1.0_dkind/(sqrt(velocity_value) + abs(div_v_int%cells(i,j,k)))	!# L2 norm

					if (delta_t_interm < time_step) then
						time_step = delta_t_interm
					end if
				end if
			end if
		end do
		end do
		end do

		!$omp end do nowait
		!$omp end parallel	
		
		
		delta_t_interm1 = time_step2
		delta_t_interm2 = time_step2
		delta_t_interm3 = time_step2

		if ((this%viscosity_flag).or.(this%diffusion_flag).or.(this%heat_trans_flag)) then
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					if (this%diffusion_flag)  then
						do spec = 1, species_number
							if (D%pr(spec)%cells(i,j,k) > 1e-10_dkind) then
								delta_t_interm1x = 1.0_dkind/4.0_dkind/(D%pr(spec)%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
								if(delta_t_interm1x < delta_t_interm1) delta_t_interm1 = delta_t_interm1x
							end if
						end do
					end if 
					if (this%viscosity_flag)	delta_t_interm2 = 1.0_dkind/4.0_dkind/(nu%cells(i,j,k)/rho%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
					if (this%heat_trans_flag)	delta_t_interm3 = 1.0_dkind/4.0_dkind/(kappa%cells(i,j,k)/1000.0_dkind/rho%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
					if (min(delta_t_interm1,delta_t_interm2) < time_step2) then
						time_step2 = min(delta_t_interm1,delta_t_interm2)
					end if
					if (min(time_step2,delta_t_interm3)	 < time_step2) then
						time_step2 = min(time_step2,delta_t_interm3)	
					end if
				end if
			end do
			end do
			end do
		end if
		
	!	this%time_step = this%initial_time_step
	!	this%time_step	=	0.25_dkind * cell_size(1) / (0.2_dkind * sqrt(10*2*Pi))
	!	this%time_step	= 1.0e-04_dkind

		print *, this%courant_fraction * time_step, time_step2
	
		this%time_step	= min(time_step2,this%courant_fraction * time_step)
		
		this%time_step	= min(5.0e-03_dkind, this%time_step)

	!	print *, this%time_step
		
		time = time + this%time_step
		
	!	this%time_step	= time_step
		
		end associate			

	end subroutine

	recursive subroutine V_cycle(this,mesh_iter,time_step,nu_1,nu_2,tol,predictor,V_cycle_depth)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: mesh_iter
		real(dkind)			,intent(in)		:: time_step
		integer				,intent(in)		:: nu_1, nu_2
		real(dkind)			,intent(in)		:: tol
		logical				,intent(in)		:: predictor
		integer				,intent(in)		:: V_cycle_depth
		
		real(dkind), dimension (3,3):: lame_coeffs
		character(len=20)			:: coordinate_system
		
		real(dkind)	,dimension(3)	:: cell_size
		integer		,dimension(3,2)	:: cons_utter_loop, cons_inner_loop, subgrid_cons_inner_loop, subgrid_cons_utter_loop, subgrid_cons_inner_loop_higher
		integer		,dimension(3)	:: offset
        real(dkind)		:: sub_F_summ, r_summ
		real(dkind)		:: beta, farfield_velocity
		
		integer	:: poisson_iteration
		logical	:: converged
		real(dkind)	:: a_norm
		
		character(len=20)		:: boundary_type_name
		
		integer	:: dimensions
		integer	:: factor
		integer	:: plus, sign, bound_number
		integer	:: i,j,k,dim,dim1,dim2
		integer	:: ilb, irb, ilt, irt, jlb, jrb, jlt, jrt
		
		

		associate (	ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					grad_F_a		=> this%grad_F_a			, &
					grad_F_b		=> this%grad_F_b			, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr		, &
					
					sub_F			=> this%sub_F				, &
					sub_R			=> this%sub_R				, &
					sub_E			=> this%sub_E				, &
					sub_E_old		=> this%sub_E_old			, &
					sub_bc			=> this%sub_bc				, &
					sub_mesh		=> this%sub_mesh)

		cons_utter_loop	= this%domain%get_local_utter_cells_bounds()	
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		
		subgrid_cons_inner_loop = this%subgrid_cons_inner_loop(mesh_iter,:,:)
		subgrid_cons_utter_loop = this%subgrid_cons_utter_loop(mesh_iter,:,:)
		if (mesh_iter < this%number_of_meshes-1) subgrid_cons_inner_loop_higher = this%subgrid_cons_inner_loop(mesh_iter+1,:,:)
		
		dimensions		= this%domain%get_domain_dimensions()
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
        
		factor		= 2**mesh_iter
		
		beta		= 2.0_dkind/3.0_dkind
		
		cell_size	= cell_size*factor
		
		converged = .false.
		
		coordinate_system	= this%domain%get_coordinate_system_name()

		if (mesh_iter < V_cycle_depth) then	
		
			!$omp parallel default(none)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs) ,			&
			!$omp& firstprivate(this) , &
			!$omp& shared(sub_E,sub_E_old,sub_R,sub_F,a_norm,beta,factor,nu_1,mesh_iter,v_f, v_f_old, F_a, F_b, coordinate_system,		&
			!$omp& cons_inner_loop,cons_utter_loop,subgrid_cons_utter_loop,subgrid_cons_inner_loop,subgrid_cons_inner_loop_higher,bc,dimensions,cell_size,mesh,converged,predictor,time_step,poisson_iteration,farfield_velocity_array,sub_F_summ,r_summ)	
					
			!$omp master
			poisson_iteration = 0
			!$omp end master
			!$omp barrier
		
			do while ((.not.converged).and.(poisson_iteration <= nu_1))
					
				!$omp barrier
				
				!$omp master
				a_norm		= 0.0_dkind
				converged	= .false.
                sub_F_summ	= 0.0_dkind
                r_summ		= 0.0_dkind
				!$omp end master
					
				!$omp barrier

				!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then

						lame_coeffs		= 1.0_dkind				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select
					
						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dkind
								
						sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) + cell_size(1)*cell_size(1)*sub_F(mesh_iter)%cells(i,j,k)								

						do dim = 1, dimensions
						
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) +	(sub_E_old(mesh_iter)%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - sub_E_old(mesh_iter)%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) -	(sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)		

							do plus = 1,2
							
								sign			= (-1)**plus
					!			bound_number	= bc%bc_markers(offset(1)+sign*I_m(dim,1),offset(2)+sign*I_m(dim,2),offset(3)+sign*I_m(dim,3))
								bound_number	= sub_bc(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
                                        !   sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))
										
                                        case('inlet')
									!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))

                                        case('outlet')	
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) - (sub_E_old(mesh_iter)%cells(i,j,k) + sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
										sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -sub_E_old(mesh_iter)%cells(i,j,k)																					
												
									end select	
								end if
							end do
						end do								
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*sub_R(mesh_iter)%cells(i,j,k)
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k))

					end if
				end do
				end do
				end do
				!$omp end do


				!$omp do collapse(3) schedule(guided)
				do k = subgrid_cons_utter_loop(3,1),subgrid_cons_utter_loop(3,2)
				do j = subgrid_cons_utter_loop(2,1),subgrid_cons_utter_loop(2,2)
				do i = subgrid_cons_utter_loop(1,1),subgrid_cons_utter_loop(1,2)
					sub_R(mesh_iter)%cells(i,j,k) =  sub_R(mesh_iter)%cells(i,j,k)/cell_size(1)/cell_size(1)
					sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
				end do
				end do
				end do			
				!$omp end do
				
				!if(a_norm/a_norm_init < tol) converged = .true.
				!$omp barrier
				!$omp master
					
				!if(a_norm/a_norm_init < tolerance) converged = .true.
				
				if ((poisson_iteration == 0).or.(poisson_iteration == nu_1)) then
				!	print *, mesh_iter, a_norm, poisson_iteration
				end if				
						
				!pause	
				
				poisson_iteration	= poisson_iteration + 1
				!$omp end master
				!$omp barrier
			
			end do
			
			!!$omp do collapse(3) schedule(guided)
			!do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
			!do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
			!do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)		
			!	offset = 1
			!	do dim = 1, dimensions
			!		offset(dim) = cons_inner_loop(dim,1) + factor * ((i-1) * I_m(dim,1) + (j-1) * I_m(dim,2) + (k-1) * I_m(dim,3))
			!	end do			
			!
			!	if(bc%bc_markers(offset(1),offset(2),offset(3)) == 0) then
			!	
			!		lame_coeffs		= 1.0_dkind				
		 !
			!		select case(coordinate_system)
			!			case ('cartesian')	
			!				lame_coeffs			= 1.0_dkind
			!			case ('cylindrical')
			!				! x -> r, y -> z
			!				lame_coeffs(1,1)	= mesh%mesh(1,offset(1),offset(2),offset(3)) - 0.5_dkind*cell_size(1)			
			!				lame_coeffs(1,2)	= mesh%mesh(1,offset(1),offset(2),offset(3))
			!				lame_coeffs(1,3)	= mesh%mesh(1,offset(1),offset(2),offset(3)) + 0.5_dkind*cell_size(1)	
			!			case ('spherical')
			!				! x -> r
			!				lame_coeffs(1,1)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)) - 0.5_dkind*cell_size(1))**2
			!				lame_coeffs(1,2)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)))**2
			!				lame_coeffs(1,3)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)) + 0.5_dkind*cell_size(1))**2
			!		end select
			!	
			!	
			!		sub_R(mesh_iter)%cells(i,j,k) = 0.0_dkind
			!					
			!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) + cell_size(1)*cell_size(1)*sub_F(mesh_iter)%cells(i,j,k)								
   !
			!		do dim = 1, dimensions
			!					
			!			sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) +	(sub_E(mesh_iter)%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - sub_E(mesh_iter)%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2) 
			!					
			!			sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) -	(sub_E(mesh_iter)%cells(i,j,k) - sub_E(mesh_iter)%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)
   !
			!			do plus = 1,2
			!				sign			= (-1)**plus
			!				bound_number	= bc%bc_markers(offset(1)+sign*I_m(dim,1),offset(2)+sign*I_m(dim,2),offset(3)+sign*I_m(dim,3))
			!				if( bound_number /= 0 ) then
			!			
			!					boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
			!					select case(boundary_type_name)
			!						case('wall','inlet')
			!							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
			!						case('outlet')	
			!							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) - sub_E(mesh_iter)%cells(i,j,k) - sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		
			!					end select		
			!				end if
			!			end do
			!		end do	
			!
			!		sub_R(mesh_iter)%cells(i,j,k) =  sub_R(mesh_iter)%cells(i,j,k)/cell_size(1)/cell_size(1)
			!	
			!	end if
   !
			!end do
			!end do
			!end do		
			!!$omp end do
			
			!$omp do collapse(3) schedule(guided) reduction(+:sub_F_summ,r_summ)
			do k = subgrid_cons_inner_loop_higher(3,1),subgrid_cons_inner_loop_higher(3,2)
			do j = subgrid_cons_inner_loop_higher(2,1),subgrid_cons_inner_loop_higher(2,2)
			do i = subgrid_cons_inner_loop_higher(1,1),subgrid_cons_inner_loop_higher(1,2)
				offset = 1
				do dim = 1, dimensions
					offset(dim) = subgrid_cons_inner_loop(dim,1) + 2 * ((i-1) * I_m(dim,1) + (j-1) * I_m(dim,2) + (k-1) * I_m(dim,3))
				end do
				
				if (dimensions == 1) then
					sub_F(mesh_iter+1)%cells(i,j,k) = 0.5_dkind*(sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)))
				end if				

				if (dimensions == 2) then
					sub_F(mesh_iter+1)%cells(i,j,k) = 0.25_dkind * (sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1),offset(2)+1,offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2)+1,offset(3)))
                end if
			
        
				select case(coordinate_system)
					case ('cartesian')	
						sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k)
                        r_summ = r_summ + 1.0_dkind! this%sub_cells_number(mesh_iter+1)
					case ('cylindrical')
						! x -> r, y -> z
						sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k) * sub_mesh(mesh_iter+1)%cells(1,i,j,k)
						r_summ		= r_summ + sub_mesh(mesh_iter+1)%cells(1,i,j,k)
                    case ('spherical')
						! x -> r
						sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k) * sub_mesh(mesh_iter+1)%cells(1,i,j,k)**2
						r_summ		= r_summ + sub_mesh(mesh_iter+1)%cells(1,i,j,k)**2
				end select

			end do
			end do
            end do
			!$omp end do
            
            if (this%all_Neumann_flag) then	
				!$omp do collapse(3) schedule(guided) 
				do k = subgrid_cons_inner_loop_higher(3,1),subgrid_cons_inner_loop_higher(3,2)
				do j = subgrid_cons_inner_loop_higher(2,1),subgrid_cons_inner_loop_higher(2,2)
				do i = subgrid_cons_inner_loop_higher(1,1),subgrid_cons_inner_loop_higher(1,2)
					sub_F(mesh_iter+1)%cells(i,j,k) = sub_F(mesh_iter+1)%cells(i,j,k) - sub_F_summ / r_summ
				end do
				end do
				end do
				!$omp end do    
            end if
            
        !    print *, sum(sub_F(mesh_iter+1)%cells(:,1,1))
            
			!$omp end parallel				
			
			call this%V_cycle(mesh_iter+1,time_step,nu_1,nu_2,tol,predictor,V_cycle_depth)

			!$omp parallel default(none)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs,ilb, irb, ilt, irt, jlb, jrb, jlt, jrt) ,			&
			!$omp& firstprivate(this) , &
			!$omp& shared(sub_E,sub_E_old,sub_R,sub_F,a_norm,beta,factor,nu_2,mesh_iter,v_f, v_f_old, F_a, F_b, coordinate_system,					&
			!$omp& cons_inner_loop,cons_utter_loop,subgrid_cons_inner_loop,subgrid_cons_utter_loop,subgrid_cons_inner_loop_higher,bc,dimensions,cell_size,mesh,converged,predictor,time_step,farfield_velocity_array,poisson_iteration)

			!$omp do collapse(3) schedule(guided)	
			do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
			do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
			do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)

				if (dimensions == 1) then
				
					ilb = int(i/2)
					irb = int(i/2)+1
					
					if (mod(i,2) == 0) then
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.25_dkind*(3.0_dkind*sub_E(mesh_iter+1)%cells(ilb,j,k) + sub_E(mesh_iter+1)%cells(irb,j,k))  
					end if			
					if (mod(i,2) /= 0) then
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.25_dkind*(sub_E(mesh_iter+1)%cells(ilb,j,k) + 3.0_dkind*sub_E(mesh_iter+1)%cells(irb,j,k))  
					end if	
				end if				
				
				if (dimensions == 2) then
				
					ilb = int(i/2)
					jlb = int(j/2)
						
					irb = int(i/2)+1
					jrb = int(j/2)
						
					ilt = int(i/2)
					jlt = int(j/2)+1
						
					irt = int(i/2)+1
					jrt = int(j/2)+1

					! Bilinear
					!if ((mod(i,2) == 0).and.(mod(j,2) /= 0)) then
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (9.0_dkind*sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(ilb,jlb,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(irt,jrt,k)) 
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (sub_E(mesh_iter+1)%cells(irb,jrb,k))
					!end if
					!if ((mod(i,2) /= 0).and.(mod(j,2) == 0)) then
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (9.0_dkind*sub_E(mesh_iter+1)%cells(irb,jrb,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(ilb,jlb,k)) 
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(irt,jrt,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					!end if
					!if ((mod(i,2) == 0).and.(mod(j,2) == 0)) then
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (9.0_dkind*sub_E(mesh_iter+1)%cells(ilb,jlb,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(irb,jrb,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (sub_E(mesh_iter+1)%cells(irt,jrt,k))
					!end if			
					!if ((mod(i,2) /= 0).and.(mod(j,2) /= 0)) then
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (9.0_dkind*sub_E(mesh_iter+1)%cells(irt,jrt,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(irb,jrb,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (3.0_dkind*sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					!	sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) + 0.0625_dkind * (sub_E(mesh_iter+1)%cells(ilb,jlb,k))
					!end if	

				! Kwak
					if ((mod(i,2) == 0).and.(mod(j,2) /= 0)) then
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k) 		+ 0.25_dkind * (2.0_dkind*sub_E(mesh_iter+1)%cells(ilt,jlt,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(ilb,jlb,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(irt,jrt,k)) 
					end if
					if ((mod(i,2) /= 0).and.(mod(j,2) == 0)) then
						sub_E_old(mesh_iter)%cells(i,j,k) =  sub_E(mesh_iter)%cells(i,j,k)		+ 0.25_dkind * (2.0_dkind*sub_E(mesh_iter+1)%cells(irb,jrb,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(ilb,jlb,k)) 
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(irt,jrt,k))
					end if
					if ((mod(i,2) == 0).and.(mod(j,2) == 0)) then
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)		+ 0.25_dkind * (2.0_dkind*sub_E(mesh_iter+1)%cells(ilb,jlb,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(irb,jrb,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					end if			
					if ((mod(i,2) /= 0).and.(mod(j,2) /= 0)) then
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)		+ 0.25_dkind *(2.0_dkind*sub_E(mesh_iter+1)%cells(irt,jrt,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(irb,jrb,k))
						sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k)	+ 0.25_dkind * (sub_E(mesh_iter+1)%cells(ilt,jlt,k))
					end if	
					
				end if				
				
			end do
			end do
			end do				
			!$omp end do
			
			!$omp master
			poisson_iteration = 0
			!$omp end master
			!$omp barrier
		
			do while ((.not.converged).and.(poisson_iteration <= nu_2))
			
				!$omp barrier
				
				!$omp master
				a_norm	= 0.0_dkind
				converged = .false.
				!$omp end master
					
				!$omp barrier
				
				!$omp do collapse(3) schedule(guided) reduction(+:a_norm)	
				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dkind				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select
					
						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dkind
								
						sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) + cell_size(1)*cell_size(1)*sub_F(mesh_iter)%cells(i,j,k)								

						do dim = 1, dimensions
								
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) +	(sub_E_old(mesh_iter)%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - sub_E_old(mesh_iter)%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2) 
								
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) -	(sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

							do plus = 1,2

								sign			= (-1)**plus
					!			bound_number	= bc%bc_markers(offset(1)+sign*I_m(dim,1),offset(2)+sign*I_m(dim,2),offset(3)+sign*I_m(dim,3))
								bound_number	= sub_bc(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
									!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))
										
                                        case('inlet')
									!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))

                                        case('outlet')	
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) - (sub_E_old(mesh_iter)%cells(i,j,k) + sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -sub_E_old(mesh_iter)%cells(i,j,k)																				
												
									end select
								end if
							end do
						end do								
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*sub_R(mesh_iter)%cells(i,j,k)
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k))
	 
					end if
				end do
				end do
				end do
				!$omp end do

				!$omp do collapse(3) schedule(guided)	
				do k = subgrid_cons_utter_loop(3,1),subgrid_cons_utter_loop(3,2)
				do j = subgrid_cons_utter_loop(2,1),subgrid_cons_utter_loop(2,2)
				do i = subgrid_cons_utter_loop(1,1),subgrid_cons_utter_loop(1,2)
					sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
				end do
				end do
				end do	
				!$omp end do
			
				!if(a_norm/a_norm_init < tol) converged = .true.
					
				!$omp barrier
				!$omp master
					
				if ((poisson_iteration == 0).or.(poisson_iteration == nu_2)) then
				!	print *, mesh_iter, a_norm, poisson_iteration
				end if	
				!pause
				!if(a_norm/a_norm_init < tolerance) converged = .true.
					
				poisson_iteration	= poisson_iteration + 1
				!$omp end master
				!$omp barrier
			
			end do
			!$omp end parallel				
	
		else
			
			!$omp parallel default(none)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs) ,			&
			!$omp& firstprivate(this) , &
			!$omp& shared(sub_E,sub_E_old,sub_R,sub_F,a_norm,beta,factor,nu_1,mesh_iter,v_f, v_f_old, F_a, F_b, coordinate_system,		&
			!$omp& cons_inner_loop,cons_utter_loop,subgrid_cons_utter_loop,subgrid_cons_inner_loop,subgrid_cons_inner_loop_higher,bc,dimensions,cell_size,mesh,converged,predictor,time_step,farfield_velocity_array,poisson_iteration)	
					
			!$omp master
			poisson_iteration = 0
			!$omp end master
			!$omp barrier
		
			do while ((poisson_iteration <= 2000*nu_1).and.(.not.converged))
					
				!$omp barrier
				
				!$omp master
				a_norm	= 0.0_dkind
				converged = .false.
				!$omp end master
					
				!$omp barrier

				!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
			
				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then

						lame_coeffs		= 1.0_dkind				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dkind
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dkind*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dkind*cell_size(1))**2
						end select

						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dkind
								
						sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) + cell_size(1)*cell_size(1)*sub_F(mesh_iter)%cells(i,j,k)								

						do dim = 1, dimensions
								
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) +	(sub_E_old(mesh_iter)%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - sub_E_old(mesh_iter)%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2) 
								
							sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) -	(sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

							do plus = 1,2
								
								sign			= (-1)**plus
					!			bound_number	= bc%bc_markers(offset(1)+sign*I_m(dim,1),offset(2)+sign*I_m(dim,2),offset(3)+sign*I_m(dim,3))
								bound_number	= sub_bc(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
									!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))
										
                                        case('inlet')
									!		sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k)	+ (sub_E_old(mesh_iter)%cells(i,j,k) - sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
                                            sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = sub_E_old(mesh_iter)%cells(i,j,k) !sub_E_old(mesh_iter)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3))

                                        case('outlet')	
											sub_R(mesh_iter)%cells(i,j,k) = sub_R(mesh_iter)%cells(i,j,k) - (sub_E_old(mesh_iter)%cells(i,j,k) + sub_E_old(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											sub_E(mesh_iter)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -sub_E_old(mesh_iter)%cells(i,j,k)																				
												
									end select		
								end if
							end do
						end do								
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*sub_R(mesh_iter)%cells(i,j,k) *beta
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k))

					end if
				end do
				end do
				end do
				!$omp end do
				
				
				!$omp do collapse(3) schedule(guided)	
				do k = subgrid_cons_utter_loop(3,1),subgrid_cons_utter_loop(3,2)
				do j = subgrid_cons_utter_loop(2,1),subgrid_cons_utter_loop(2,2)
				do i = subgrid_cons_utter_loop(1,1),subgrid_cons_utter_loop(1,2)
					sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
				end do
				end do
				end do			
				!$omp end do
				

				
				!$omp barrier
				!$omp master
					
				!	print *, a_norm, poisson_iteration
				!	pause
				!print *, sqrt(a_norm), poisson_iteration
				!pause
				!if(a_norm/a_norm_init < tolerance) converged = .true.
					
				if(a_norm < 1.0e-05) then
                    converged = .true.
                    print *, poisson_iteration
                end if
				poisson_iteration	= poisson_iteration + 1

				!$omp end master
				!$omp barrier
			
			end do
		
		!	print *, mesh_iter, a_norm, poisson_iteration
			
		end if
		
		!$omp end parallel

		continue
		
		end associate
		
	end subroutine V_cycle
	
	subroutine farfield_values_modifier(this, time)
		class(fds_solver)	,intent(in)		:: this
		real(dkind)			,intent(in)		:: time
		
		integer				:: j, bound
		character(len=20)	:: boundary_type_name
		real(dkind)			:: RND 
		real(dkind)	,save	:: farfield_velocity
		integer		,save	:: iteration
		
		iteration = iteration + 1 
		
		if (iteration == 1) then
			do bound = 1, size(this%boundary%bc_ptr%boundary_types)
				boundary_type_name = this%boundary%bc_ptr%boundary_types(bound)%get_type_name()
				if (boundary_type_name == 'inlet') then
					farfield_velocity = this%boundary%bc_ptr%boundary_types(bound)%get_farfield_velocity()
				end if
			end do
		end if

		!if (mod(iteration,100) == 0) then
		!	do j = 1,size(farfield_velocity_array)
		!		call RANDOM_NUMBER(RND)  
		!		farfield_velocity_array(j) = farfield_velocity + RND*(0.5_dkind*farfield_velocity)
		!	end do
		!end if
		
		if (time < 5e-03) then
			farfield_velocity_array = 0.0_dkind
		else
			farfield_velocity_array = farfield_velocity
		end if	
	end subroutine
	
	subroutine igniter(this, time)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time
		
		real(dkind)				:: duration, delay
		integer					:: i,j,k
		integer	,dimension(3,2)	:: cons_inner_loop
		integer		,save		:: iteration = 0

		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		!delay	= 5.0e-03_dkind
		delay	= 0.0e-03_dkind
		duration = 5.0e-03_dkind
		
		associate (	rho	=> this%rho%s_ptr)
			if ((time > delay).and. (iteration == 0)) then !(time <= delay + duration)) then !(time <= delay + relax)) then
				iteration = iteration + 1
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if ((i-128.5)**2 + (j-40)**2 < 157) then
						rho%cells(i,j,:)	= this%rho_0 - this%rho_0*4.0_dkind/5.0! *(time-delay)/relax
						!T%cells(i,j,:)	= 1500.0_dkind
					end if
				end do
				end do
				end do
				
				call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.true.)	
			end if
			
		!	if (time > delay + relax) pause
		end associate
		
	end subroutine	
	
	subroutine write_data_table(this,table_file)
		class(fds_solver)	,intent(inout)  :: this	
		character(len=*)	,intent(in)		:: table_file
		
		real(dkind)			:: max_T, min_T, max_grad_T, lf
		integer				:: max_grad_index
		
		character(len=100)						:: name_string
		character(len=10)						:: fmt
		real(dkind)	,dimension(:), allocatable	:: temp

		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		real	,dimension(3)	:: cell_size
		
		integer	:: io_unit
		integer :: table_size
		integer	:: i, spec, dim
		integer	:: ierr
		
		open(newunit = io_unit	, file = trim(task_setup_folder) // trim(fold_sep) // trim(table_file), status = 'replace'	, form = 'formatted') 	
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		allocate(temp(species_number+dimensions+3))
		
		associate (	rho				=> this%rho%s_ptr			, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

		max_grad_T = 0.0	
		max_T = 0.0
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1) > max_grad_T) then
				max_grad_T = abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1)
				max_grad_index = i
			end if
			if(T%cells(i,1,1) > max_T)	max_T = T%cells(i,1,1)
		end do

		min_T = 300.0
		
		lf = (max_T-min_T)/max_grad_T
		
		name_string = 'x'
		name_string = trim(name_string) // '  ' // trim(T%name_short)
		name_string = trim(name_string) // '  ' // trim(rho%name_short)

		do dim = 1,dimensions
			name_string = trim(name_string) // '  ' // trim(v%pr(dim)%name_short)
		end do			
		
		do spec = 1,species_number
			name_string = trim(name_string) // '  ' // trim(Y%pr(spec)%name_short)
		end do		
		
		write(fmt,'(A,I2,A)') '(', species_number+dimensions+3, 'E14.7)'
		
		write(io_unit,'(A)',iostat = ierr)	name_string
		
		do i = max_grad_index-10*lf/cell_size(1), max_grad_index+10*lf/cell_size(1)
			temp(1) = mesh%mesh(1,i,1,1)
			temp(2) = T%cells(i,1,1)
			temp(3) = rho%cells(i,1,1)
			do dim = 1,dimensions
				temp(dim+3) = v%pr(dim)%cells(i,1,1)
			end do
			do spec = 1,species_number
				temp(spec+dimensions+3) = Y%pr(spec)%cells(i,1,1)
			end do			
			write(io_unit,fmt,iostat = ierr) temp
		end do	
			
		print *, max_grad_index, lf 
			
		end associate
		
	end subroutine	
	
	
	pure function get_time_step(this)
		real(dkind)						:: get_time_step
		class(fds_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dkind)						:: get_time
		class(fds_solver)	,intent(in)		:: this

		get_time = this%time
	end function

	
	
	
end module
	
!subroutine Weighted_Jacobi_solver_MG(cell_size,max_iter,u,r,f,bc,boundary_types)
!
!	use kind_parameters
!	use global_data
!	use boundary_type_class
!
!    real(dkind)								,intent(in)     :: cell_size
!    real(dkind)			,dimension(:,:,:)	,intent(inout)  :: u, r, f, bc
!	type(boundary_type)	,dimension(:)		,intent(in)		:: boundary_types
!    integer									,intent(in)     :: max_iter
!
!    real(dkind)			,dimension(:,:,:)	,allocatable    :: u_new
!
!	integer	,dimension(3,2)	:: bound_sl
!    integer			:: cells_number(3), dimensions
!    integer			:: i,j,k, dim,iter
!    real(dkind)	    :: tolerance, omega, residual_sum, residual_sum_init
!
!	tolerance = 1e-02_dkind
!    omega = 2.0_dkind/3.0_dkind
!
!	dimensions = 0
!	do dim = 1, 3
!		cells_number(dim) = size(sub_u(mesh)%cells,dim)
!		if (cells_number(dim) /= 1) dimensions = dimensions + 1
!	end do
!	
!	do dim = 1, 3
!		bound_sl(dim,1) = lbound(sub_u(mesh)%cells,dim)
!		bound_sl(dim,2) = ubound(sub_u(mesh)%cells,dim)
!	end do		
!	
!    allocate(u_new(cells_number(1),cells_number(2),cells_number(3)))
!    u_new = u
!	
!	do k = bound_sl(3,1),bound_sl(3,2)
!	do j = bound_sl(2,1),bound_sl(2,2)
!	do i = bound_sl(1,1),bound_sl(1,2)
!		r(i,j,k) = f(i,j,k)
!		do dim = 1, dimensions
!			if( bc(i,j,k) == 0) then
!				r(i,j,k) = r(i,j,k) - (u(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))-2.0_dkind*u(i,j,k)+u(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/(cell_size**2)
!			else
!				bound_number	= bc(i,j,k)
!				if( bound_number /= 0 ) then
!					boundary_type_name = boundary_types(bound_number)%get_type_name()
!					select case(boundary_type_name)
!						case('wall')
!							r(i,j,k) = 
!						case('outlet')
!							r(i,j,k) = 
!					end select
!				end if
!			end if
!		end do
!    end do
!	end do
!	end do
!
!    residual_sum_init = sum(r)
!    residual_sum      = residual_sum_init
!
!    iter = 0
!    do while ((iter < max_iter).and.(residual_sum/residual_sum_init > tolerance))
!		do k = bound_sl(3,1),bound_sl(3,2)
!		do j = bound_sl(2,1),bound_sl(2,2)
!		do i = bound_sl(1,1),bound_sl(1,2)
!			if( bc(i,j,k) == 0) then
!				u_new(i,j,k) = - (cell_size**2 * f(i,j,k))/(2.0_dkind*dimensions)
!			
!				do dim = 1, dimensions
!					u_new(i,j,k) =  u_new(i,j,k) + (u(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + u(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/(2.0_dkind*dimensions) 
!				end do
!            
!				u_new(i,j,k) = (1.0 - omega)*u(i,j,k) + omega*u_new(i,j,k)
!			else
!				bound_number	= bc(i,j,k)
!				if( bound_number /= 0 ) then
!					boundary_type_name = boundary_types(bound_number)%get_type_name()
!					select case(boundary_type_name)
!						case('wall')
!						case('outlet')
!					end select
!				end if
!			end if
!		end do
!		end do
!		end do
!
!        do i = 2,N-1
!            r(i) = f(i) - (u_new(i-1) - 2.0*u_new(i) + u_new(i+1))/(cell_size**2)
!        end do
!
!        u = u_new
!
!        residual_sum = sum(r)
!        iter = iter + 1
!    end do
!
!end subroutine Weighted_Jacobi_solver_MG
!
!recursive subroutine V_cycle(mesh, sub_u, sub_r, sub_f, nu_1, nu_2, cell_size)
!
!    use subgrid_type
!
!    implicit none
!
!    type(subgrid)   ,dimension(:)  ,intent(inout)   :: sub_u, sub_r, sub_f
!    integer                        ,intent(in)      :: nu_1, nu_2
!    integer                        ,intent(in)      :: mesh
!    real                           ,intent(in)      :: cell_size
!
!    integer     :: number_of_meshes
!    integer     :: cells_number(3), offset
!    integer     :: i, j, k, dim
!
!    number_of_meshes = size(sub_u)
!
!	do dim = 1, 3
!		cells_number(dim) = size(sub_u(mesh)%cells,dim)
!	end do
!	
!    call Weighted_Jacobi_solver_MG(cell_size, nu_1, sub_u(mesh)%cells, sub_r(mesh)%cells, sub_f(mesh)%cells)
!
!    if (mesh /= 1) then
!
!		do j = 1, 2**(mesh-1)+1
!		do k = 1, 2**(mesh-1)+1
!        do i = 1, 2**(mesh-1)+1
!			!    offset = 1 + 2**(number_of_meshes - (mesh-1))*(i-1)
!			offset_x = 1 + 2*(i-1)  
!			offset_y = 1 + 2*(j-1)
!			offset_z = 1 + 2*(k-1)
!			sub_f(mesh-1)%cells(i,j,k) = sub_r(mesh)%cells(offset_x, offset_y, offset_z)
!		end do
!		end do
!		end do
!
!        call V_cycle(mesh-1, sub_u, sub_r, sub_f, nu_1, nu_2, cell_size*2)
!
!        do i = 1, cells_number
!		do j = 1, cells_number
!		do k = 1, cells_number
!            if (mod(i,2) == 0).and.(mod(j,2) /= 0) then
!                sub_u(mesh)%cells(i,j,k) = sub_u(mesh)%cells(i,j,k) + 0.5*(sub_u(mesh-1)%cells(i/2,int(j/2),k) + sub_u(mesh-1)%cells(i/2+1,int(j/2),k))  
!            end if
!            if (mod(i,2) /= 0).and.(mod(j,2) == 0) then
!                sub_u(mesh)%cells(i,j,k) = sub_u(mesh)%cells(i,j,k) + 0.5*(sub_u(mesh-1)%cells(int(i/2),j/2,k) + sub_u(mesh-1)%cells(int(i/2),j/2+1,k))  
!            end if
!            if (mod(i,2) == 0).and.(mod(j,2) == 0) then
!                sub_u(mesh)%cells(i,j,k) = sub_u(mesh)%cells(i,j,k) + 0.25*(sub_u(mesh-1)%cells(i/2,j/2,k) + sub_u(mesh-1)%cells(i/2+1,j/2,k) + sub_u(mesh-1)%cells(i/2,j/2+1,k) + sub_u(mesh-1)%cells(i/2+1,j/2+1,k))  
!            end if			
!			if (mod(i,2) /= 0).and.(mod(j,2) /= 0) then
!                sub_u(mesh)%cells(i) = sub_u(mesh)%cells(i) + sub_u(mesh-1)%cells((i+1)/2,(j+1)/2)
!            end if
!        end do    
!		end do
!		end do
!    end if
!
!    call Weighted_Jacobi_solver_MG(cell_size, nu_1, sub_u(mesh)%cells, sub_r(mesh)%cells, sub_f(mesh)%cells)
!
!end subroutine
