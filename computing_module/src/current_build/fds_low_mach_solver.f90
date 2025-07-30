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
    
    use benchmarking
    use supplementary_routines
	
	implicit none
	
#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: fds_solver, fds_solver_c

	type(field_scalar_cons)	,target	:: p_dyn	,p_stat		,p_stat_old	,dp_stat_dt	,p_int, T_int, rho_int, rho_old	
	type(field_scalar_cons)	,target	:: div_v	,div_v_int	,ddiv_v_dt	,H	, H_old, R
	type(field_scalar_cons)	,target	:: E_f_int
	type(field_vector_cons)	,target	:: v_int, Y_int, Y_old, v_prod_sources
	type(field_scalar_flow)	,target	:: F_a, F_b
	type(field_vector_flow)	,target	:: v_f, v_f_old	

	real(dp)	,dimension(:)	,allocatable	:: concs
	real(dp)	,dimension(:)	,allocatable	:: Y_rho_int
	!$omp threadprivate(concs,Y_rho_int)

    type(timer)     :: fds_timer
    type(timer)     :: fds_gas_dynamics_timer
    type(timer)     :: fds_multigrid_timer
    type(timer)     :: fds_eos_timer
    type(timer)     :: fds_chemistry_timer
    type(timer)     :: fds_diffusion_timer
    type(timer)     :: fds_heattransfer_timer
    type(timer)     :: fds_viscosity_timer
  
	real(dp)	,dimension(:)	,allocatable	:: farfield_velocity_array
	integer	:: flame_structure_unit

	type subgrid 
		real(dp), dimension(:,:,:), allocatable :: cells
    end type
    
	type subgrid_vector
		real(dp), dimension(:,:,:,:), allocatable :: cells
	end type

	type fds_solver
		logical			:: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, hydrodynamics_flag, CFL_condition_flag, all_Neumann_flag, perturbed_velocity
		real(dp)		:: courant_fraction
		real(dp)		:: time, time_step, initial_time_step
        real(dp)    , dimension(3)  :: g
		integer			:: additional_particles_phases_number, additional_droplets_phases_number
		
		logical			:: perturbed_velocity_flag, stabilizing_inlet_flag, igniter_flag 
        
		type(viscosity_solver)				:: visc_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(diffusion_solver)				:: diff_solver
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(table_approximated_real_gas)	:: state_eq

		type(lagrangian_droplets_solver), dimension(:)	    ,allocatable	:: droplets_solver			!# Lagrangian droplets solver
!		type(droplets_solver)			, dimension(:)	    ,allocatable	:: droplets_solver			!# Continuum droplets solver
		
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
		
		type(field_vector_cons_pointer)	:: v		, v_prod_gd		, v_prod_visc	, v_prod_sources	, v_int	
		type(field_vector_cons_pointer)	:: Y		, Y_prod_diff	, Y_prod_chem	, Y_int				, Y_old
		type(field_vector_cons_pointer)	:: D
		type(field_vector_flow_pointer)	:: v_f		, v_f_old
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	::  rho_prod_droplets, E_f_prod_droplets, E_f_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  Y_prod_droplets		
		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Lagrangian droplets solver
!		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets				!# Continuum droplets solver
		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Lagrangian particles solver
!		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  v_prod_particles			!# Continuum particles solver		


		real(dp)	,dimension(:,:,:,:)	,allocatable	:: vorticity, grad_F_a, grad_F_b
		real(dp)	,dimension(:,:,:)	,allocatable	:: p_old
		
		integer												:: number_of_meshes
		integer			,dimension(:,:,:)	,allocatable	:: subgrid_cons_inner_loop, subgrid_cons_utter_loop
		type(subgrid)	,dimension(:)		,allocatable	:: sub_E, sub_F, sub_R, sub_E_old, sub_bc
        type(subgrid_vector)	,dimension(:)		,allocatable	:: sub_mesh
        real(dp)				,dimension(:)		,allocatable	:: sub_cells_number
        real(dp)				,dimension(:)		,allocatable	:: sub_cell_size
		
        integer					:: load_counter
        
		real(dp)				:: rho_0
		real(dp)             	:: calc_time
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

		real(dp)						:: calculation_time		
		
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
        
		real(dp)	,dimension(3)	:: cell_size
        integer     ,dimension(3)   :: offset
		real(dp)					:: x, y
		real(dp)					:: farfield_velocity
		integer	:: dimensions
		integer	:: i,j,k,n,dim,dim1
        
        integer								    		:: neighbours
        integer,        dimension(:,:),	    allocatable	:: neighbours_indexes
        integer										    :: resid
        logical                                         :: inside_flag

        integer	:: counter
		integer	:: mesh
		
		integer :: cells_number
		
		constructor%calc_time           = 0.0_dp
		constructor%diffusion_flag		= problem_solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= problem_solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= problem_solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= problem_solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= problem_solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= problem_solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= problem_solver_options%get_CFL_condition_flag()
        
        !# sub solver options
		    constructor%perturbed_velocity_flag	= .false.
            constructor%stabilizing_inlet_flag	= .false.
		    constructor%igniter_flag	        = .false.
            
        constructor%g                       = problem_solver_options%get_grav_acc()
        
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

        if(constructor%perturbed_velocity) then
			call manager%create_vector_field(v_prod_sources	,'velocity_production_sources'		,'v_prod_sources'	,'spatial')
			constructor%v_prod_sources%v_ptr	=> v_prod_sources
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
				constructor%droplets_solver(droplets_phase_counter)	= lagrangian_droplets_solver_c(manager, droplets_params, droplets_phase_counter)		!# Lagrangian droplets solver
!				constructor%droplets_solver(droplets_phase_counter)	= droplets_solver_c(manager, droplets_params, droplets_phase_counter)					!# Continuum droplets solver
				write(var_name,'(A,I2.2)') 'energy_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_droplets', droplets_phase_counter						
				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)								!# Lagrangian droplets solver
				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr						!# Lagrangian droplets solver
!				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)						!# Continuum droplets solver								
!				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Continuum droplets solver
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

		constructor%vorticity= 0.0_dp
		constructor%p_old	 = 0.0_dp
		constructor%grad_F_a = 0.0_dp									
		constructor%grad_F_b = 0.0_dp

		cons_utter_loop	= manager%domain%get_local_utter_cells_bounds()	
		cons_inner_loop = manager%domain%get_local_inner_cells_bounds()	
		flow_inner_loop	= manager%domain%get_local_inner_faces_bounds()	
		
		problem_data_io				= data_io_c(manager,calculation_time)									
		
        if(problem_data_io%get_load_counter() /= 0) then
			call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
			call problem_data_io%add_io_scalar_cons_field(constructor%div_v)
		end if
        
		call problem_data_io%input_all_data()

        constructor%load_counter	= problem_data_io%get_load_counter()
        
		dimensions		= manager%domain%get_domain_dimensions()

		do dim = 1, dimensions		
			loop = flow_inner_loop
			do dim1 = 1,dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do
			do k = loop(3,1),loop(3,2)
			do j = loop(2,1),loop(2,2) 
			do i = loop(1,1),loop(1,2)
				do dim1 = 1, dimensions
					constructor%v_f%v_ptr%pr(dim1)%cells(dim,i,j,k) = 0.5_dp * (constructor%v%v_ptr%pr(dim1)%cells(i,j,k) + constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) )
				end do			
			end do
			end do
			end do
        end do   
        
        constructor%p_stat%s_ptr%cells	= constructor%p%s_ptr%cells        
        
		if(constructor%load_counter == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
			call problem_data_io%add_io_scalar_cons_field(constructor%div_v)
		
			constructor%p_dyn%s_ptr%cells		= 0.0_dp 
	
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
        
 		if(constructor%load_counter == 1) then
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
        else
			call constructor%state_eq%apply_state_equation_low_mach_fds(problem_solver_options%get_initial_time_step(),predictor=.true.)
			call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
        end if	       

        cell_size						= constructor%mesh%mesh_ptr%get_cell_edges_length()
		
		constructor%time				= calculation_time
		constructor%initial_time_step	= problem_solver_options%get_initial_time_step()
		constructor%time_step			= constructor%initial_time_step

		constructor%rho_0				= constructor%rho%s_ptr%cells(cons_inner_loop(1,2),cons_inner_loop(2,2) ,1)!constructor%rho%s_ptr%cells(1,1,1)
!		print *,constructor%rho_0

		species_number = manager%chemistry%chem_ptr%species_number
		
		do bound = 1, size(constructor%boundary%bc_ptr%boundary_types)
			boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound)%get_type_name()
			if (boundary_type_name == 'inlet') then
				farfield_velocity = constructor%boundary%bc_ptr%boundary_types(bound)%get_farfield_velocity()
			end if
		end do		

		allocate(farfield_velocity_array(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))

        if(constructor%load_counter == 1) then
			farfield_velocity_array = farfield_velocity
        else
            farfield_velocity_array = constructor%v%v_ptr%pr(1)%cells(1,cons_inner_loop(2,2)/2,1)
		end if
            
		!$omp parallel
		allocate(concs(species_number))
		allocate(Y_rho_int(species_number))
		concs		= 0.0_dp
		Y_rho_int	= 0.0_dp
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
			
														

			constructor%sub_E(mesh)%cells		= 0.0_dp
			constructor%sub_F(mesh)%cells		= 0.0_dp
			constructor%sub_R(mesh)%cells		= 0.0_dp
			constructor%sub_E_old(mesh)%cells	= 0.0_dp
			constructor%sub_bc(mesh)%cells		= 0.0_dp
            constructor%sub_mesh(mesh)%cells	= 0.0_dp
            
            do dim = 1, dimensions
				cells_number = cells_number * (constructor%subgrid_cons_inner_loop(mesh,dim,2) - constructor%subgrid_cons_inner_loop(mesh,dim,1) + 1)
            end do
		
            constructor%sub_cells_number(mesh) = cells_number
            constructor%sub_cell_size(mesh) = cell_size(1)*2**mesh
		end do	
		
		constructor%sub_bc(0)%cells	= constructor%boundary%bc_ptr%bc_markers
        constructor%sub_mesh(0)%cells	=  constructor%mesh%mesh_ptr%mesh
		
		constructor%all_Neumann_flag = .true.
		
        do k = constructor%subgrid_cons_utter_loop(0,3,1),constructor%subgrid_cons_utter_loop(0,3,2)
		do j = constructor%subgrid_cons_utter_loop(0,2,1),constructor%subgrid_cons_utter_loop(0,2,2)
		do i = constructor%subgrid_cons_utter_loop(0,1,1),constructor%subgrid_cons_utter_loop(0,1,2)
			bound_number = constructor%boundary%bc_ptr%bc_markers(i,j,k)
            if (bound_number /= 0) then
			    boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound_number)%get_type_name()
									
			    if ((boundary_type_name == 'outlet')) then
				    constructor%all_Neumann_flag = .false.
                end if
            end if
		end do
		end do
        end do
        
		do mesh = 0,constructor%number_of_meshes-2
	
            do k = constructor%subgrid_cons_inner_loop(mesh+1,3,1),constructor%subgrid_cons_inner_loop(mesh+1,3,2)
			do j = constructor%subgrid_cons_inner_loop(mesh+1,2,1),constructor%subgrid_cons_inner_loop(mesh+1,2,2)
			do i = constructor%subgrid_cons_inner_loop(mesh+1,1,1),constructor%subgrid_cons_inner_loop(mesh+1,1,2)            
				offset = 1
				do dim = 1, dimensions
					offset(dim) = constructor%subgrid_cons_inner_loop(mesh,dim,1) + 2 * ((i-1) * I_m(dim,1) + (j-1) * I_m(dim,2) + (k-1) * I_m(dim,3))
				end do
				
				do dim = 1, dimensions
					constructor%sub_mesh(mesh+1)%cells(dim,i,j,k) = 0.5_dp*(constructor%sub_mesh(mesh)%cells(dim,offset(1),offset(2),offset(3)) + constructor%sub_mesh(mesh)%cells(dim,offset(1)+1,offset(2),offset(3)))
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
		
        
        neighbours = 2**dimensions
        allocate(neighbours_indexes(neighbours,3))
        
        
        do mesh = 1,constructor%number_of_meshes-1

			do k = constructor%subgrid_cons_utter_loop(mesh,3,1),constructor%subgrid_cons_utter_loop(mesh,3,2)
			do j = constructor%subgrid_cons_utter_loop(mesh,2,1),constructor%subgrid_cons_utter_loop(mesh,2,2)
			do i = constructor%subgrid_cons_utter_loop(mesh,1,1),constructor%subgrid_cons_utter_loop(mesh,1,2)
                
                neighbours_indexes	= 1
                do n = 1, neighbours
                    resid = n - 1
					do dim = dimensions, 1, -1
						neighbours_indexes(n,dim) = -int(resid / 2**(dim-1)) + int(2*i)*I_m(dim,1) + int(2*j)*I_m(dim,2) + int(2*k)*I_m(dim,3)
                        resid = mod(resid, 2**(dim-1))
                    end do
                end do

                counter = 0
                do n = 1, neighbours
                    inside_flag = .true.
                    do dim = 1, dimensions
                        if ((constructor%subgrid_cons_utter_loop(mesh-1,dim,1) <= neighbours_indexes(n,dim)).and.(neighbours_indexes(n,dim) <= constructor%subgrid_cons_utter_loop(mesh-1,dim,2))) then
                            inside_flag = inside_flag .and. .true.
                        else
                            inside_flag = inside_flag .and. .false.
                        end if
                    end do
                    
                    if (inside_flag) then
                        counter = counter + 1
                        constructor%sub_bc(mesh)%cells(i,j,k) = constructor%sub_bc(mesh)%cells(i,j,k) + constructor%sub_bc(mesh-1)%cells(neighbours_indexes(n,1),neighbours_indexes(n,2),neighbours_indexes(n,3))
					end if
                end do
                
                constructor%sub_bc(mesh)%cells(i,j,k) = int(constructor%sub_bc(mesh)%cells(i,j,k)/(counter))

			end do
			end do
			end do
            
        end do
        
        
        !!!!!######
		!constructor%sub_bc(5)%cells(2,0,1) = 3.0

        call RANDOM_SEED()
        
		continue
		
    call manager%create_timer(fds_timer                 ,'FDS solver time'              , 'sol_t')
    call manager%create_timer(fds_gas_dynamics_timer    ,'FDS gas dynamics time'        , 'gd_t')
    call manager%create_timer(fds_eos_timer             ,'FDS eos solver time'          , 'eos_t')
    call manager%create_timer(fds_chemistry_timer       ,'FDS chemistry solver time'    , 'chem_t')
    call manager%create_timer(fds_diffusion_timer       ,'FDS diffusion solver time'    , 'diff_t')
    call manager%create_timer(fds_heattransfer_timer    ,'FDS heattransfer solver time' , 'ht_t')
    call manager%create_timer(fds_viscosity_timer       ,'FDS viscosity solver time'    , 'visc_t')
    call manager%create_timer(fds_multigrid_timer       ,'Multigrid solver time'        , 'mg_t')
 
	end function

	subroutine solve_problem(this,iteration,stop_flag)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: iteration
        logical				,intent(inout)	:: stop_flag
		
		integer	:: droplets_phase_counter, particles_phase_counter
		integer	:: specie

		logical	:: stabilized_flag

        call fds_timer%tic()
        
		this%time = this%time + this%time_step		

		if (this%igniter_flag) then
			call this%igniter(this%time)
        end if
		
!		call this%farfield_values_modifier(this%time)
		
!		if (iteration > 10) perturbed_velocity_field = .false.

        call fds_chemistry_timer%tic()
		if (this%reactive_flag)				call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
        call fds_chemistry_timer%toc(new_iter=.true.)
        
        call fds_diffusion_timer%tic()
		if (this%diffusion_flag)			call this%diff_solver%solve_diffusion(this%time_step)
        call fds_diffusion_timer%toc(new_iter=.true.)
        
        call fds_viscosity_timer%tic()
		if (this%viscosity_flag)			call this%visc_solver%solve_viscosity(this%time_step)
        call fds_viscosity_timer%toc(new_iter=.true.)
        
        call fds_heattransfer_timer%tic()
		if (this%heat_trans_flag)			call this%heat_trans_solver%solve_heat_transfer(this%time_step)				
        call fds_heattransfer_timer%toc(new_iter=.true.)
 
		if (this%perturbed_velocity_flag)		call this%perturb_velocity_field(this%time_step)
        
		if(this%additional_droplets_phases_number /= 0) then
			do droplets_phase_counter = 1, this%additional_droplets_phases_number
				call this%droplets_solver(droplets_phase_counter)%droplets_solve(this%time_step)				!# Lagrangian droplets solver
				!call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_main(this%time)    !# Continuum droplets solver
				!call this%droplets_solver(droplets_phase_counter)%droplets_euler_step_v_E(this%time_step)		!# Continuum droplets solver
				!call this%droplets_solver(droplets_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum droplets solver
				!call this%droplets_solver(droplets_phase_counter)%droplets_lagrange_step(this%time_step)		!# Continuum droplets solver
				!call this%droplets_solver(droplets_phase_counter)%droplets_final_step(this%time_step)			!# Continuum droplets solver		
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
		
        call fds_gas_dynamics_timer%tic()
		call this%calculate_interm_Y_predictor(this%time_step)
        call fds_gas_dynamics_timer%toc()
        
        call fds_eos_timer%tic()
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.true.)		
        call fds_eos_timer%toc(new_iter=.true.)
        
        call fds_gas_dynamics_timer%tic()
		call this%apply_boundary_conditions(this%time_step,predictor=.true.)		
		call this%calculate_divergence_v		(this%time_step,predictor=.true.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.true.)
		call this%calculate_velocity			(this%time_step,predictor=.true.)
        call fds_gas_dynamics_timer%toc(new_iter=.true.)
		        
        call fds_viscosity_timer%tic()
		if (this%viscosity_flag)		call this%visc_solver%solve_viscosity(this%time_step)
        call fds_viscosity_timer%toc()
		if (this%perturbed_velocity_flag)	call this%perturb_velocity_field(this%time_step)
        
        call fds_gas_dynamics_timer%tic()
		call this%calculate_interm_Y_corrector(this%time_step)
		call fds_gas_dynamics_timer%toc()
		
        call fds_eos_timer%tic()
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.false.)
		call fds_eos_timer%toc(new_iter=.true.)		
        
        call fds_gas_dynamics_timer%tic()
		call this%apply_boundary_conditions(this%time_step,predictor=.false.)
		call this%calculate_divergence_v		(this%time_step,predictor=.false.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.false.)
		call this%calculate_velocity			(this%time_step,predictor=.false.)		

        if (this%CFL_condition_flag) then
			call this%calculate_time_step()
		end if
        
        call fds_gas_dynamics_timer%toc(new_iter=.true.)
 
!		call this%chem_kin_solver%write_chemical_kinetics_table('15_pcnt_H2-Air_table(T).dat')

		if (this%stabilizing_inlet_flag) then
			call this%stabilizing_inlet_1D(this%time, stabilized_flag)
			if (stabilized_flag) then
				stop_flag = .true.
!				call this%chem_kin_solver%write_chemical_kinetics_table('15.0_pcnt_H2-Air_table(T).dat')
!				call this%write_data_table('H2-Air_flamelet.dat')
				stop
			end if
		end if
         
        
        call fds_timer%toc(new_iter=.true.)
 
		!call this%state_eq%check_conservation_laws()

	end subroutine

	subroutine calculate_interm_Y_predictor(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		
		real(dp)	:: B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dp)	:: spec_summ
		
		real(dp)	,dimension(3)	:: cell_size		
		real(dp)	:: rhs
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		real(dp), dimension (3,3)	:: lame_coeffs
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
                    rho_old         => this%rho_old%s_ptr		, &
					v_f				=> this%v_f%v_ptr			, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_droplets	=> this%Y_prod_droplets	, &
                    Y_old           => this%Y_old%v_ptr         , &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)
					
			!$omp parallel default(shared)  private(i,j,k,dim,spec,spec_summ,rhs,flux_right,flux_left,lame_coeffs) !, &
			!!$omp& firstprivate(this)
            !!$omp& shared(cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,time_step)					

            !$omp do collapse(3) schedule(static)	
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	

					lame_coeffs		= 1.0_dp				
				
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dp
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
					end select					
				
					do spec = 1, species_number
					
						Y_rho_int(spec) = rho%cells(i,j,k)*Y%pr(spec)%cells(i,j,k) 
						
						rhs = 0.0_dp
						
						do dim = 1, dimensions
													
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
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
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
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
			
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					spec_summ = 0.0_dp
					do spec = 1,species_number
						spec_summ = spec_summ + max(Y_int%pr(spec)%cells(i,j,k), 0.0_dp)
					end do
					do spec = 1,species_number
						Y_int%pr(spec)%cells(i,j,k) = max(Y_int%pr(spec)%cells(i,j,k), 0.0_dp) / spec_summ
					end do
				end if
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(static)
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
            end associate
			
		!	print *, '1',Y_int%pr(1)%cells(50,20,1),rho_int%cells(50,20,1)
			
			continue
			


	end subroutine
	
	subroutine calculate_divergence_v(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
	
		real(dp)	:: flux_left, flux_right, specie_enthalpy, mixture_cp, dp_dt, D_sum, P_sum, U_sum, div_sum, average_molar_mass, mol_mix_conc
		
		real(dp)					:: cell_volume
		real(dp)	,dimension(3)	:: cell_size, cell_surface_area
			
		real(dp)					:: energy_source = 1.0e06_dp
		real(dp)	,save			:: time
		integer		,save			:: iter = 1

		real(dp), dimension (3,3)	:: lame_coeffs
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
		
        D_sum = 0.0_dp
        P_sum = 0.0_dp		
        U_sum = 0.0_dp
        
        div_sum = 0.0_dp
        
        time = time + 0.5_dp*time_step
        
        if(time <= 1.0e-02_dp) then
            energy_source = 2.0e02_dp
        else
            energy_source = 0.0_dp
        end if
        
    !	print *, '2', div_v_int%cells(50,20,1),T%cells(50,20,1),D_sum,P_sum
		
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

			!$omp parallel default(shared)  private(flux_right,flux_left,i,j,k,dim,spec,mixture_cp,specie_enthalpy,plus,sign,bound_number,cell_volume,cell_surface_area,lame_coeffs, average_molar_mass, mol_mix_conc) !, &
			!!$omp& firstprivate(this)
            !!$omp& shared(D_sum,P_sum,U_sum,energy_source,cons_inner_loop,species_number,dimensions,cell_size,coordinate_system)					
 
            !$omp do collapse(3) schedule(static)	reduction(+:D_sum,P_sum)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
				
					cell_volume	= mesh%get_cell_volume()
					lame_coeffs	= 1.0_dp

					select case(coordinate_system)
						case ('cartesian')	

						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
							cell_volume			= cell_volume * mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
							cell_volume			= cell_volume * mesh%mesh(1,i,j,k)**2
					end select					

					div_v_int%cells(i,j,k) = 0.0_dp
					
					!if ((i > cons_inner_loop(1,1) + 1).and.(j > cons_inner_loop(2,1) + 1).and.(i < cons_inner_loop(1,2) - 1).and.(j < cons_inner_loop(2,2) - 1)) then
					if (this%viscosity_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_visc%cells(i,j,k)	![J/m^3/s]
					if (this%heat_trans_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_heat%cells(i,j,k)	![J/m^3/s]
					if (this%diffusion_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_diff%cells(i,j,k)	![J/m^3/s]
					if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_chem%cells(i,j,k)	![J/m^3/s]
					!end if

					average_molar_mass = 0.0_dp
					do spec = 1,species_number
						average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
					end do
				
					mol_mix_conc		= 1.0_dp / average_molar_mass

					concs = 0.0_dp
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
							if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
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
							if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
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
					
						specie_enthalpy = (thermo%calculate_specie_enthalpy(T%cells(i,j,k),spec) - thermo%calculate_specie_enthalpy(T_ref,spec))  / thermo%molar_masses(spec)
					
						do dim = 1, dimensions
					
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
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
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
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
					P_sum = P_sum + (1.0_dp / p_stat%cells(i,j,k) - 1.0_dp / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc)) * cell_volume 
				end if
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(static)	reduction(+:U_sum)
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
										if(dim==2) cell_surface_area(dim) = cell_surface_area(dim) * mesh%mesh(1,i,j,k)			! - 0.5_dp*cell_size(1)
									case ('spherical')
										! x -> r
										if(dim==1) cell_surface_area(dim) = cell_surface_area(dim) * (mesh%mesh(1,i,j,k))**2		! - 0.5_dp*cell_size(1)
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

			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)			
				if(bc%bc_markers(i,j,k) == 0) then	
					
					average_molar_mass = 0.0_dp
					do spec = 1,species_number
						average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
                    end do

					mol_mix_conc	= 1.0_dp / average_molar_mass				
				
                    concs = 0.0_dp
					do spec = 1,species_number
						concs(spec)				= Y%pr(spec)%cells(i,j,k) *  mol_mix_conc / thermo%molar_masses(spec)
					end do	    
      
					mixture_cp		= thermo%calculate_mixture_cp(T%cells(i,j,k), concs)
				
					dp_stat_dt%cells(i,j,k) = (D_sum - U_sum)/ P_sum

					if (this%all_Neumann_flag) then	
						div_v_int%cells(i,j,k)  = div_v_int%cells(i,j,k) - (1.0_dp / p_stat%cells(i,j,k) - 1.0_dp / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mol_mix_conc))* dp_stat_dt%cells(i,j,k)
					end if

				end if
			end do
			end do
			end do			
			!$omp end do

			!$omp end parallel

            end associate
		
			iter = iter + 1
		!	print *, '2', div_v_int%cells(50,20,1),T%cells(50,20,1),D_sum,P_sum
			
		continue
	end subroutine	
	
	subroutine calculate_pressure_poisson(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		
		real(dp)	:: H_center, H_left, H_right, F_a_left, F_a_right, F_b_left, F_b_right, sub_F_summ, r_summ1, r_summ2
		real(dp)	:: H_residual, H_max, H_max_old, H_average, residual, a_norm_init, a_norm, a_norm_prev, a_normfinal, H_summ, sum_ddiv_v_dt, grad_F_a_summ, grad_F_b_summ
		real(dp)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dp), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		integer		:: r_i, r_j, r_k

		real(dp)	,dimension(3)	:: cell_size		
		real(dp)	:: time_step_adj, beta, beta_old, spectral_radii_jacobi, edges_number
		real(dp), save	:: time = 0.0_dp
		
		integer	:: dimensions, iterations, iterations_number
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop, subgrid_cons_inner_loop
		integer	,dimension(3,2)	:: loop
		integer	,dimension(3)	:: offset
		
		character(len=20)		:: boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3
		real(dp)				:: bc_coeff
		
		integer	:: droplets_phase_counter, particles_phase_counter
		integer	:: sign, bound_number, bound_number1, bound_number2, bound_number3
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		integer			:: poisson_iteration, pressure_iteration, overall_poisson_iteration, mesh_iter, v_cycle_iteration
		integer			:: nu_0, nu_1, nu_2, V_cycle_depth
		real(dp)		:: tolerance
		
		logical			:: converged = .false., v_cycle_converged = .false.
		logical			:: pressure_converged = .false.
		
		dimensions		= this%domain%get_domain_dimensions()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		subgrid_cons_inner_loop	= this%subgrid_cons_inner_loop(1,:,:)
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		time_step_adj	=  0.2_dp*cell_size(1)**2								!	0.025_dp*(3.0_dp * cell_size(1)**2) !0.1_dp*cell_size(1)**2
		
		if (predictor) time = time + time_step

            associate (	vorticity		=> this%vorticity			, &
                        ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
					    F_a				=> this%F_a%s_ptr			, &
					    F_b				=> this%F_b%s_ptr)
                
			ddiv_v_dt%cells	= 0.0_dp
            F_a%cells		= 0.0_dp		
			F_b%cells		= 0.0_dp	
            vorticity		= 0.0_dp
            
            end associate
            
			sum_ddiv_v_dt	= 0.0_dp
			grad_F_a_summ	= 0.0_dp
			r_summ1			= 0.0_dp
            r_summ2			= 0.0_dp
            
            associate ( rho_old			=> this%rho_old%s_ptr		, &
                        rho_int			=> this%rho_int%s_ptr		, &
                        div_v_int		=> this%div_v_int%s_ptr		, &
                        ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                        F_a				=> this%F_a%s_ptr			, &
                        vorticity		=> this%vorticity			, &
                        grad_F_a		=> this%grad_F_a			, &
                        v_f				=> this%v_f%v_ptr			, &
                        v_f_old			=> this%v_f_old%v_ptr		, &
                        v_prod_visc		=> this%v_prod_visc%v_ptr	, &
                        v_prod_droplets	=> this%v_prod_droplets		, &
                        v_prod_particles=> this%v_prod_particles	, &
                        v_prod_sources	=> this%v_prod_sources%v_ptr , &
                        mesh			=> this%mesh%mesh_ptr		, &
                        bc				=> this%boundary%bc_ptr		)

			!$omp parallel default(shared)  private(i,j,k,dim,dim1,dim2,loop,lame_coeffs,farfield_velocity,sign,bound_number,bound_number1,bound_number2,bound_number3,plus,boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3,bc_coeff) !, &
			!!$omp& firstprivate(this)
            !!$omp& shared(sum_ddiv_v_dt,grad_F_a_summ,cons_inner_loop,cons_utter_loop,flow_inner_loop,dimensions,predictor,cell_size,coordinate_system,time_step,r_summ)					
                        
            !$omp do collapse(3) schedule(static)	reduction(+:sum_ddiv_v_dt,r_summ1) 
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
				
					lame_coeffs		= 1.0_dp				
		
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dp
                            
                            r_summ1				= r_summ1 + 1.0_dp
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
                            
                            r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                            
                            r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)**2
					end select					
				
					if (predictor) then
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / time_step
						
						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k) - (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3)  - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2) / time_step
						end do
					else
						ddiv_v_dt%cells(i,j,k)	= div_v_int%cells(i,j,k) / (0.5_dp * time_step)
						
						do dim = 1, dimensions
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k)  - 0.5_dp * ( (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)	&
																							+ (v_f_old%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f_old%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)) / (0.5_dp * time_step)
						end do
					end if

					sum_ddiv_v_dt = sum_ddiv_v_dt + ddiv_v_dt%cells(i,j,k) * lame_coeffs(1,2)
				end if
			end do
			end do
			end do	
			!$omp end do

			if (this%all_Neumann_flag) then
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					ddiv_v_dt%cells(i,j,k)	=	ddiv_v_dt%cells(i,j,k) - sum_ddiv_v_dt/r_summ1
				end if
			end do
			end do
			end do		
			!$omp end do			
			end if
			
			!$omp do collapse(3) schedule(static)
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
					if(abs(vorticity(dim,i,j,k)) < 1e-03) vorticity(dim,i,j,k) = 0.0_dp
				end do
				
				
			end do
			end do
			end do
			!$omp end do
			
			!$omp do collapse(3) schedule(static)	
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
											!	vorticity(3,i,j,k) = vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) !0.0_dp
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = vorticity(3,i,j,k) !0.0_dp
											!end if
										case('inlet')
											!if (sign < 0) then
											!	vorticity(3,i,j,k) = 0.0_dp
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp
											!end if	
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i,j,k) = 0.0_dp !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)) !0.5_dp*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))) !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) =  cell_size(1)*(vorticity(3,i,j,k) - (v_f%pr(dim)%cells(dim,i,j,k) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) = -cell_size(1)*(vorticity(3,i,j,k) + (v_f%pr(dim)%cells(dim,i,j,k) - v_f%pr(dim)%cells(dim,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i,j,k)/cell_size(1)) 
												end if
											else
												vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) !0.5_dp*(vorticity(3,i,j,k) + vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
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
											!	vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) = -vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) !0.0_dp
											!else
											!	vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = -vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) !0.0_dp
											!end if
										case('inlet')
											!if (sign < 0) then
											!	vorticity(3,i,j,k) = 0.0_dp
											!else
											!	vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dp
											!end if	
											continue
										case('outlet')
											if (sign < 0) then
												vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) = 0.0_dp !2.0*vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) - vorticity(3,i+2*I_m(dim,1)+I_m(dim2,1),j+2*I_m(dim,2)+I_m(dim2,2),k+2*I_m(dim,3)+I_m(dim2,3)) !vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))!2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3)) !0.5_dp*(vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))) !2.0*vorticity(3,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - vorticity(3,i+2*I_m(dim,1),j+2*I_m(dim,2),k+2*I_m(dim,3))
												if (dim2 == 1) then
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) =  cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1))
												else
													v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) = -cell_size(1)*(vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) + (v_f%pr(dim)%cells(dim,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - v_f%pr(dim)%cells(dim,i,j,k))/cell_size(1) - v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))/cell_size(1)) 
												end if
											else
												vorticity(3,i+I_m(dim,1)+I_m(dim2,1),j+I_m(dim,2)+I_m(dim2,2),k+I_m(dim,3)+I_m(dim2,3)) = 0.0_dp !2.0*vorticity(3,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3)) - vorticity(3,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) !vorticity(3,i,j,k)!2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) !0.5_dp*(vorticity(3,i,j,k) + vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) !2.0*vorticity(3,i,j,k) - vorticity(3,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
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
				
				!$omp do collapse(3) schedule(static)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then 
						if (predictor) then
							if(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) > 1e-010) then
								F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* this%g(dim))
							end if

							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(0.5_dp*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
							if (this%perturbed_velocity)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_sources%pr(dim)%cells(i,j,k) + v_prod_sources%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
                        else
							if(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) > 1e-010) then
								F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) *(this%rho_0 - rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* this%g(dim))
							end if
								
							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) *(0.5_dp*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
							if (this%perturbed_velocity)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_sources%pr(dim)%cells(i,j,k) + v_prod_sources%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
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
				
				!$omp do collapse(3) schedule(static)
				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)
				
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then 
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
                                F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)	![m/s^2]						!# Lagrangian droplets solver
								if (bc%bc_markers(i,j,k) == 0) then
!									F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum droplets solver
								end if
							end do		
						end if
							
						if (this%additional_particles_phases_number /= 0) then
							do particles_phase_counter = 1, this%additional_particles_phases_number
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)	![m/s^2]						!# Lagrangian particles solver
								!if (bc%bc_markers(i,j,k) == 0) then
								!	F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum particles solver
								!end if
							end do		
						end if							
							
						do dim1 = 1, 3
						do dim2 = 1, dimensions
								if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
								if(((dim-dim1) == 1).or.((dim-dim1) == -2)) then
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - 0.25_dp * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								else
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.25_dp * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								end if
								end if
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

!$omp do collapse(3) schedule(static) reduction(+:grad_F_a_summ, r_summ2)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					
					lame_coeffs		= 1.0_dp				
			
					select case(coordinate_system)
						case ('cartesian')	
							lame_coeffs			= 1.0_dp
                            
                            r_summ2				= r_summ2 + 1.0_dp
						case ('cylindrical')
							! x -> r, y -> z
							lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
							lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
							lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
                            
                            r_summ2				= r_summ2 + mesh%mesh(1,i,j,k)
						case ('spherical')
							! x -> r
							lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
							lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
							lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                            
                            r_summ2				= r_summ2 + mesh%mesh(1,i,j,k)**2
					end select						
					
					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)	= (F_a%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_a%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
						grad_F_a_summ		= grad_F_a_summ + grad_F_a(dim,i,j,k) * lame_coeffs(1,2)
					end do

				end if
			end do
			end do
			end do		
			!$omp end do			
			
!            !$omp barrier
!            print *, 'R summ', r_summ 
            
			if (this%all_Neumann_flag) then	
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
						grad_F_a(dim,i,j,k)  = grad_F_a(dim,i,j,k)  - grad_F_a_summ/r_summ2
					end do
				end if
			end do
			end do
			end do		
			!$omp end do			
			end if

			!$omp end parallel
			
            end associate
			     
			overall_poisson_iteration = 0			
			pressure_iteration	= 0
			pressure_converged	= .false.
			
			do while ((.not.pressure_converged).and.(pressure_iteration < 200)) 

                associate (	p_dyn				=> this%p_dyn%s_ptr , &
                            p_old				=> this%p_old)
                    
                    p_old	= p_dyn%cells
                    
                end associate				

				H_max_old	= 10.0
				a_norm_init = 0.0_dp

				grad_F_b_summ = 0.0_dp
                
                r_summ1 = 0.0_dp
				
				associate ( p_dyn			=> this%p_dyn%s_ptr			, &
                            rho_old			=> this%rho_old%s_ptr		, &
                            rho_int			=> this%rho_int%s_ptr		, &
                            ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                            F_a				=> this%F_a%s_ptr			, &
                            F_b				=> this%F_b%s_ptr			, &
                            H_old			=> this%H_old%s_ptr			, &
                            R				=> this%R%s_ptr				, &

                            grad_F_a		=> this%grad_F_a			, &
                            grad_F_b		=> this%grad_F_b			, &
                            v				=> this%v%v_ptr				, &
                            v_f				=> this%v_f%v_ptr			, &
                            v_f_old			=> this%v_f_old%v_ptr		, &

                            mesh			=> this%mesh%mesh_ptr		, &
                            bc				=> this%boundary%bc_ptr		)

				!$omp parallel default(shared)  private(i,j,k,dim,dim2,loop,plus,sign,bound_number,boundary_type_name,residual,lame_coeffs,farfield_velocity) !, &
                !!$omp& firstprivate(this)
                !!$omp& shared(grad_F_b_summ,predictor,cons_utter_loop,cons_inner_loop,dimensions,cell_size,coordinate_system,a_norm_init,time_step,farfield_velocity_array,r_summ)                

				do dim = 1, dimensions
					loop(3,1) = cons_inner_loop(3,1)
					loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

					loop(2,1) = cons_inner_loop(2,1)
					loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

					loop(1,1) = cons_inner_loop(1,1)
					loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))	

					!$omp do collapse(3) schedule(static)
					do k = loop(3,1),loop(3,2)
					do j = loop(2,1),loop(2,2)
					do i = loop(1,1),loop(1,2)

						if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then
							if (predictor) then	
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_old%cells(i,j,k))		&
														/	(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_old%cells(i,j,k))	&
														*	(1.0_dp/rho_old%cells(i,j,k)	- 1.0_dp/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	!/ cell_size(1)
							else
								F_b%cells(dim,i,j,k)=	-	(p_dyn%cells(i,j,k)	*rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_int%cells(i,j,k))		&
														/	(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_int%cells(i,j,k))	&
														*	(1.0_dp/rho_int%cells(i,j,k)	- 1.0_dp/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	!/ cell_size(1)
							end if

						end if
					end do
					end do
					end do
					!$omp end do

				end do
				
				!$omp do collapse(3) schedule(static) reduction(+:grad_F_b_summ,r_summ1)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dp				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
                                
                                r_summ1				= r_summ1 + 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
                                
                                r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
                                
                                r_summ1				= r_summ1 + mesh%mesh(1,i,j,k)**2
						end select						
					
						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)	= (F_b%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_b%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
							grad_F_b_summ	= grad_F_b_summ + grad_F_b(dim,i,j,k) * lame_coeffs(1,2)
						end do
						
					end if
				end do
				end do
				end do		
				!$omp end do			
			
				if (this%all_Neumann_flag) then	
				!$omp do collapse(3) schedule(static)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1, dimensions
							grad_F_b(dim,i,j,k)  = grad_F_b(dim,i,j,k)  - grad_F_b_summ/r_summ1
						end do
					end if
				end do
				end do
				end do		
				!$omp end do				
				end if
				
				!$omp do collapse(3) schedule(static) reduction(+:a_norm_init)
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dp				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select						
					
						R%cells(i,j,k)	= 0.0_dp
						
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
											R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																			- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																					+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
										
											!R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
                                        case('inlet')
										!	farfield_velocity = farfield_velocity_array(factor * j)
											farfield_velocity = farfield_velocity_array(1)
											if(predictor) then
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																				- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																							
											else
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																				- sign*(	farfield_velocity - 0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																						+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
											end if
										case('outlet')
											R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
											if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dp) then
												do dim2 = 1, dimensions
													R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
												end do
											end if
									end select
								end if
							end do
						end do

						a_norm_init = a_norm_init + abs(R%cells(i,j,k)*lame_coeffs(1,2))
						
					end if
				end do
				end do
				end do
				!$omp end do

				!$omp end parallel	
				
                end associate

				poisson_iteration	= 0 
				converged			= .false.
				H_max				= 0.0_dp
				a_norm_prev			= a_norm_init
				beta				= 2.0_dp/3.0_dp
				v_cycle_converged	= .false.
			
				nu_0 = 5
				nu_1 = 2
				nu_2 = 2
				V_cycle_depth = this%number_of_meshes - 2

				tolerance = 1e-02_dp
					
                call fds_multigrid_timer%tic()
					
				v_cycle_iteration = 0
				do while (((v_cycle_iteration <= 0).or.(.not.v_cycle_converged)).and.(v_cycle_iteration <= nu_0))
			
					nu_1 = 2!(v_cycle_iteration+2) * 2
					nu_2 = 1!(v_cycle_iteration+2) * 2
				
					associate (     ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                                    F_a				=> this%F_a%s_ptr			, &
                                    F_b				=> this%F_b%s_ptr			, &
                                    H				=> this%H%s_ptr				, &
                                    H_old			=> this%H_old%s_ptr			, &
                                    R				=> this%R%s_ptr				, &
                                    grad_F_a		=> this%grad_F_a			, &
                                    grad_F_b		=> this%grad_F_b			, &
                                    v				=> this%v%v_ptr				, &
                                    v_f				=> this%v_f%v_ptr			, &
                                    v_f_old			=> this%v_f_old%v_ptr		, &
                                    mesh			=> this%mesh%mesh_ptr		, &
                                    bc				=> this%boundary%bc_ptr		, &
                        
                                    sub_F			=> this%sub_F				, &
                                    sub_E			=> this%sub_E)
					
					do while (poisson_iteration <= nu_1)
					
						a_norm	= 0.0_dp
						converged = .false.
					
                        !$omp parallel default(shared)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) !, &
                        !!$omp& firstprivate(this)
				                        
						!$omp do collapse(3) schedule(static) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
			  
								lame_coeffs		= 1.0_dp				
		
								select case(coordinate_system)
									case ('cartesian')	
										lame_coeffs			= 1.0_dp
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
								end select								

								R%cells(i,j,k) = 0.0_dp
								
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
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
													!R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
												
                                                case('inlet')
												!	farfield_velocity = farfield_velocity_array(factor * j)
													farfield_velocity = farfield_velocity_array(1)
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																						- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																							
													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																						- sign*(	farfield_velocity - 0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																								+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
													R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
													if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dp) then
														do dim2 = 1, dimensions
															R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
														end do
													end if										
											end select		
										end if
									end do
								end do								
								
								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*R%cells(i,j,k)
  
								a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))
							end if
						end do
						end do
						end do
						!$omp end do
					
						!$omp do collapse(3) schedule(static) 
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							R%cells(i,j,k)	= R%cells(i,j,k) / cell_size(1) / cell_size(1)
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do
                        !$omp end parallel
						
      !                  if ((poisson_iteration == 0).or.(poisson_iteration == nu_1)) then
      !                      if (poisson_iteration == 0) write(*,'(A,/,A,/,A)') '*****', 'Starting Multigrid', '*****'
						!	print *, "Initial V cycle step", a_norm, poisson_iteration
						!end if

						poisson_iteration	= poisson_iteration + 1
					end do

					a_norm		= 0.0_dp
					converged   = .false.
                    sub_F_summ	= 0.0_dp
                    r_summ1		= 0.0_dp
				
                    !$omp parallel default(shared)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) !, &
                    !!$omp& firstprivate(this)

					!$omp do collapse(3) schedule(static) reduction(+:sub_F_summ,r_summ1)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                        if(bc%bc_markers(i,j,k) == 0) then
							sub_F(0)%cells(i,j,k)	= R%cells(i,j,k)
                        
							select case(coordinate_system)
								case ('cartesian')	
									sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k)
									r_summ1 = r_summ1 + 1.0_dp! this%sub_cells_number(mesh_iter+1)
								case ('cylindrical')
									! x -> r, y -> z
									sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k) * mesh%mesh(1,i,j,k)
									r_summ1		= r_summ1 + mesh%mesh(1,i,j,k)
								case ('spherical')
									! x -> r
									sub_F_summ	= sub_F_summ + sub_F(0)%cells(i,j,k) * mesh%mesh(1,i,j,k)**2
									r_summ1		= r_summ1 + mesh%mesh(1,i,j,k)**2
							end select
                        end if    
                      
					end do
					end do
                    end do
					!$omp end do					

                    if (this%all_Neumann_flag) then	
                        !$omp do collapse(3) schedule(static) 
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                        if(bc%bc_markers(i,j,k) == 0) then
							sub_F(0)%cells(i,j,k) = sub_F(0)%cells(i,j,k) - sub_F_summ / r_summ1 !this%sub_cells_number(0)
                        end if
					end do
					end do
                    end do
					!$omp end do    
                    end if
                    
					!$omp end parallel
                    end associate
					
					call this%V_cycle(0,time_step,nu_1,nu_2,tolerance,predictor,V_cycle_depth)
		
					associate (     ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                                    F_a				=> this%F_a%s_ptr			, &
                                    F_b				=> this%F_b%s_ptr			, &
                                    H				=> this%H%s_ptr				, &
                                    H_old			=> this%H_old%s_ptr			, &
                                    R				=> this%R%s_ptr				, &
                                    grad_F_a		=> this%grad_F_a			, &
                                    grad_F_b		=> this%grad_F_b			, &
                                    v				=> this%v%v_ptr				, &
                                    v_f				=> this%v_f%v_ptr			, &
                                    v_f_old			=> this%v_f_old%v_ptr		, &
                                    mesh			=> this%mesh%mesh_ptr		, &
                                    bc				=> this%boundary%bc_ptr		, &

                                    sub_E			=> this%sub_E				)

					!$omp parallel default(shared)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) !,			&
					!!$omp& firstprivate(this)

					!$omp do collapse(3) schedule(static)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						H_old%cells(i,j,k) = H%cells(i,j,k) + sub_E(0)%cells(i,j,k)
					end do
					end do
					end do
					!$omp end do					

					!$omp do collapse(3) schedule(static) reduction(+:a_norm)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						if(bc%bc_markers(i,j,k) == 0) then
			  
							lame_coeffs		= 1.0_dp				
		
							select case(coordinate_system)
								case ('cartesian')	
									lame_coeffs			= 1.0_dp
								case ('cylindrical')
									! x -> r, y -> z
									lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
									lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
									lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
								case ('spherical')
									! x -> r
									lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
									lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
									lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
							end select								

							R%cells(i,j,k) = 0.0_dp
								
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
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
												
                                                !R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
											
                                            case('inlet')
											!	farfield_velocity = farfield_velocity_array(factor * j)
												farfield_velocity = farfield_velocity_array(1)
												if(predictor) then
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																					- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
												else
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																					- sign*(	farfield_velocity - 0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																							+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
												end if
											case('outlet')
												R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
												if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dp) then
													do dim2 = 1, dimensions
														R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
													end do
												end if
										end select		
									end if
								end do
							end do	
							
							a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))
							
						end if
					end do
					end do
					end do			
					!$omp end do					
                    !$omp end parallel	
					
					poisson_iteration = 0
					
					do while (poisson_iteration <= nu_2)
					
						a_norm	= 0.0_dp
						converged = .false.

						!$omp parallel default(shared)  private(i,j,k,dim,dim2,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity,offset) !,			&
					    !!$omp& firstprivate(this)

						!$omp do collapse(3) schedule(static) reduction(+:a_norm)
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
			  
								lame_coeffs		= 1.0_dp				
		
								select case(coordinate_system)
									case ('cartesian')	
										lame_coeffs			= 1.0_dp
									case ('cylindrical')
										! x -> r, y -> z
										lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
										lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
										lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
									case ('spherical')
										! x -> r
										lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
										lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
										lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
								end select								

								R%cells(i,j,k) = 0.0_dp
								
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
													R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																					- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																							+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
													!R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
												
                                                case('inlet')
												!	farfield_velocity = farfield_velocity_array(factor * j)
													farfield_velocity = farfield_velocity_array(1)
													if(predictor) then
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																						- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
													else
														R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		&
																						- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))											&
																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)				&
																						- sign*(	farfield_velocity - 0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))	&
																								+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if
												case('outlet')
													R%cells(i,j,k) = R%cells(i,j,k) - H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) 
													if ( sign*v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) >= 0.0_dp) then
														do dim2 = 1, dimensions
															R%cells(i,j,k) = R%cells(i,j,k) + v%pr(dim2)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))**2
														end do
													end if
											end select		
										end if
									end do
								end do
								
								H%cells(i,j,k)	= H_old%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*R%cells(i,j,k)
  
								a_norm = a_norm + abs(R%cells(i,j,k)*lame_coeffs(1,2))
	 
								R%cells(i,j,k)	= R%cells(i,j,k) / cell_size(1) / cell_size(1)
								
							end if
						end do
						end do
						end do
						!$omp end do
					
						!$omp do collapse(3) schedule(static) 
						do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
						do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
						do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
							H_old%cells(i,j,k) = H%cells(i,j,k)
						end do
						end do
						end do
						!$omp end do
						
						!if(a_norm/a_norm_init < tolerance) converged = .true.
					
					!$omp end parallel						
					
					    if ((poisson_iteration == 0).or.(poisson_iteration == nu_2)) then
							print *, "Final V cycle step", a_norm, poisson_iteration
                            if (poisson_iteration == nu_2) write(*,'(A,/,A,/,A)') '*****', 'Finalizing Multigrid', '*****'
						end if

                        poisson_iteration	= poisson_iteration + 1

                    end do				
                    end associate
                    
                    v_cycle_iteration = v_cycle_iteration + 1
				
					do mesh_iter = 0, this%number_of_meshes-1
						this%sub_E_old(mesh_iter)%cells	= 0.0_dp
						this%sub_E(mesh_iter)%cells		= 0.0_dp
					end do	
					
					if(abs(a_norm_init) > 1e-10_dp) then
					!	print *,a_norm,a_norm_init,a_norm/a_norm_init
						if((a_norm/a_norm_init < tolerance).or.(a_norm < 100e-0)) then
							v_cycle_converged = .true.
						end if	
					end if
				end do
				
                call fds_multigrid_timer%toc(new_iter=.true.)
				
				print *, 'Poisson iteration:', poisson_iteration

				overall_poisson_iteration = overall_poisson_iteration + poisson_iteration

				associate (	    p_dyn			=> this%p_dyn%s_ptr			, &
					            F_a				=> this%F_a%s_ptr			, &
					            F_b				=> this%F_b%s_ptr			, &
					            H				=> this%H%s_ptr				, &
					            v				=> this%v%v_ptr				, &
					            v_f				=> this%v_f%v_ptr			, &
					            v_f_old			=> this%v_f_old%v_ptr		, &
					            mesh			=> this%mesh%mesh_ptr		, &
					            bc				=> this%boundary%bc_ptr		)
                
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
											H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  H%cells(i,j,k) &
																											- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																													+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)		
										case('outlet')
											if (sign == 1) then		!# Правая граница
												if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dp) then		!# Выток, берутся значения слева
												
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
							
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/farfield_density + 0.5_dp*(farfield_velocity **2)
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	- H%cells(i,j,k) 
													!do dim2 = 1, dimensions
													!	H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + v%pr(dim2)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))**2
													!end do
													
												end if
											else					!# Левая граница
												if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dp) then										!# Выток, берутся значения справа
										 
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
												
													H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/farfield_density + 0.5_dp*(farfield_velocity **2) 
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
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)	&	
																												- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)/time_step
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i,j,k)														&
																												- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																														+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)	&	
																												- sign*(	farfield_velocity - 0.5_dp*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																														+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))) * 1.0_dp * cell_size(1)/(0.5_dp*time_step)
											end if
									end select
								end if
							end do	
						end do
					end if
				end do
				end do
				end do	
				
                end associate
                
				call this%calculate_dynamic_pressure(time_step,predictor)
				pressure_converged = .true.				
				
                associate (     p_old           => this%p_old               , &
                                p_dyn			=> this%p_dyn%s_ptr			, &
                                rho_old			=> this%rho_old%s_ptr		, &
                                rho_int			=> this%rho_int%s_ptr		, &
                                bc				=> this%boundary%bc_ptr		)

				!$omp parallel default(shared)  private(i,j,k,dim,H_residual,r_i,r_j,r_k,plus,sign,bound_number,boundary_type_name,farfield_pressure) !, &
				!!$omp& firstprivate(this)
                !!$omp& shared(H_summ,predictor,pressure_converged,cons_inner_loop,cons_utter_loop,dimensions,cell_size,time_step,time)

                !$omp do collapse(3) schedule(static)	reduction(.and.:pressure_converged)	
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
					if(bc%bc_markers(i,j,k) == 0) then
						do dim = 1,dimensions
							if (predictor) then
								if (pressure_converged) then
								if (abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dp/rho_old%cells(i,j,k)	- 1.0_dp/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>20.0_dp/cell_size(1)/cell_size(1)) then
									pressure_converged = .false.
									print *, 'Pressure error', abs(((p_dyn%cells(i,j,k) - p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dp/rho_old%cells(i,j,k)	- 1.0_dp/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))
								end if
								end if
							else
								if (pressure_converged) then
								if (abs(((p_dyn%cells(i,j,k) - p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dp/rho_int%cells(i,j,k)	- 1.0_dp/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>20.0_dp/cell_size(1)/cell_size(1)) then
									pressure_converged = .false.
									print *, 'Pressure error',  abs(((p_dyn%cells(i,j,k) - p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dp/rho_int%cells(i,j,k)	- 1.0_dp/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))
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

                end associate

				pressure_iteration = pressure_iteration + 1

				continue
			end do
				
			print *, 'Pressure iteration:', pressure_iteration
			print *, 'Overall pressure iteration:', overall_poisson_iteration
            print *, ' '
			
			continue
			
	end subroutine
	
	subroutine calculate_dynamic_pressure(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor	
				
		real(dp)	:: farfield_density, farfield_pressure, farfield_velocity
		real(dp)	:: vel_abs
		
		real(dp)	,dimension(3)	:: cell_size		
		
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
                    p_dyn			=> this%p_dyn%s_ptr			, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					H				=> this%H%s_ptr				, &
					v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					bc				=> this%boundary%bc_ptr)
	
			!$omp parallel default(shared)  private(i,j,k,dim,vel_abs,plus,sign,bound_number,boundary_type_name,farfield_density,farfield_pressure) !, &
			!!$omp& firstprivate(this)
            !!$omp& shared(predictor,cons_inner_loop,dimensions,cell_size)	

			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then

					vel_abs				= 0.0_dp
					do dim = 1,dimensions
						vel_abs				= vel_abs + (0.5_dp*(v_f%pr(dim)%cells(dim,i,j,k) + v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))))**2	
					end do	
					
					if (predictor) then	
						p_dyn%cells(i,j,k)	=  (H%cells(i,j,k) - 0.5_dp * vel_abs)*rho_old%cells(i,j,k)
					else
						p_dyn%cells(i,j,k)	=  (H%cells(i,j,k) - 0.5_dp * vel_abs)*rho_int%cells(i,j,k)
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

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*(	v%pr(dim)%cells(i,j,k) **2)  
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if
                                    case('outlet')
										!if (sign == 1) then		!# Правая граница
										!	if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dp) then		!# Выток, берутся значения слева
										!		p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	else																					!# Вток, берутся значения справа
										!		p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	end if
										!else					!# Левая граница
										!	if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dp) then										!# Выток, берутся значения справа
										!		p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	else																					!# Вток, берутся значения справа
										!		p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i,j,k)
										!	end if
										!end if
											
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*(	v%pr(dim)%cells(i,j,k) **2)   
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if                                        
                                        
									case('inlet')
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*(	v%pr(dim)%cells(i,j,k) **2)   
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
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor

		real(dp)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions, iterations
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		integer					:: plus, sign, bound_number
		character(len=20)		:: boundary_type_name
        
        real(dp)				:: farfield_velocity
       
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

		!$omp parallel default(shared)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,loop) !, &

			do dim = 1, dimensions
				loop(3,1) = cons_inner_loop(3,1)
				loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))

				loop(2,1) = cons_inner_loop(2,1)
				loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	

				loop(1,1) = cons_inner_loop(1,1)
				loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))		
				
                !$omp do collapse(3) schedule(static)
				do k = loop(3,1), loop(3,2)		
				do j = loop(2,1), loop(2,2)		
				do i = loop(1,1), loop(1,2)		
				                    
					if((bc%bc_markers(i,j,k) == 0).or.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

						if (predictor)	then
							v_f_old%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k)
							v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) - time_step * (F_a%cells(dim,i,j,k) + F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(1))
						else
							v_f%pr(dim)%cells(dim,i,j,k) = 0.5_dp * (v_f%pr(dim)%cells(dim,i,j,k) + v_f_old%pr(dim)%cells(dim,i,j,k)) - (0.5_dp * time_step) * (F_a%cells(dim,i,j,k) + F_b%cells(dim,i,j,k)  +  (H%cells(i,j,k) - H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/cell_size(1) )
                        end if	
                        
                    end if
				end do
				end do
				end do
                !$omp end do
 
			end do
				
            !$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)	
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1, dimensions
                        v%pr(dim)%cells(i,j,k) = 0.5_dp * ( v_f%pr(dim)%cells(dim,i,j,k) + v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
					end do	
				end if
			end do
			end do
            end do
			!$omp end do
                    
            !$omp do collapse(3) schedule(static)
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
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			= -v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											end if
										end do
									case('outlet')
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  v%pr(dim1)%cells(i,j,k)
											end if
										end do
									case('inlet')
										farfield_velocity		=  farfield_velocity_array(j)

										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  farfield_velocity
											else
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  0.0_dp 
											end if
										end do
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
            
		!	print *, '4', v_f%pr(2)%cells(2,50,20,1), v_f%pr(1)%cells(1,50,20,1)
			
			continue
		end associate
	end subroutine
	
	subroutine perturb_velocity_field(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step

		real(dp)	,dimension(3)	:: cell_size		
		real(dp) ,save			:: time = 0.0_dp
        
		integer	:: dimensions, iterations
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus

		real(dp)	:: kappa_turb, gamma_turb, kx_turb, ky_turb, kappa_turb_max, time_scale
		
		dimensions		= this%domain%get_domain_dimensions()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
		flow_inner_loop	= this%domain%get_local_inner_faces_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		associate (	v_f				=> this%v_f%v_ptr			, &
					v_prod_sources	=> this%v_prod_sources%v_ptr	, &
					bc				=> this%boundary%bc_ptr)

            time_scale		= 1.0e-03_dp
            kappa_turb_max	= 7.0e05_dp
            
            kappa_turb	= min(kappa_turb_max, kappa_turb_max/time_scale * time)
			if (time > 50.0e-03) kappa_turb = 0.0_dp
            
			call RANDOM_NUMBER(gamma_turb)
			
			gamma_turb = 2.0_dp * pi * gamma_turb   !2.0_dp*gamma_turb - 1.0_dp   !cos(alpha) sin(alpha)=sqrt(1-gamma_turb**2.0)
			
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                do dim = 1, dimensions
					if((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then	

						!call RANDOM_NUMBER(gamma_turb)
						!gamma_turb = 2.0_dp*gamma_turb - 1.0_dp
						!v_f%pr(dim)%cells(dim,i,j,k) = v_f%pr(dim)%cells(dim,i,j,k) + kappa_turb*gamma_turb*sqrt(this%time_step)
							
                        !gamma_turb  = atan(i*cell_size(1)/j*cell_size(1))
                        
						kx_turb = 2.0_dp*pi*sin(gamma_turb)/0.001_dp   !0.002   !sqrt(1.0_dp-gamma_turb**2.0)/0.002  
						ky_turb = 2.0_dp*pi*cos(gamma_turb)/0.001_dp   !0.002   !gamma_turb/0.002
							
						if(dim == 1)then
							v_prod_sources%pr(dim)%cells(i,j,k) = kappa_turb*cos(gamma_turb)*sqrt(time_step)*cos(kx_turb*(i-0.5)*cell_size(1)+ky_turb*(j-0.5)*cell_size(1))
						endif
						if(dim == 2)then
							v_prod_sources%pr(dim)%cells(i,j,k) = -kappa_turb*sin(gamma_turb)*sqrt(time_step)*cos(kx_turb*(i-0.5)*cell_size(1)+ky_turb*(j-0.5)*cell_size(1))
						endif                        
					end if
				end do
			end do
			end do
            end do

            time = time + time_step
            
			continue
		end associate
	end subroutine
	
	
	subroutine calculate_interm_Y_corrector(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		
		real(dp)	:: rhs, B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dp)	:: spec_summ
		
		real(dp), dimension (3,3):: lame_coeffs
		
		real(dp)	,dimension(3)	:: cell_size		
		
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
            		rho_old			=> this%rho_old%s_ptr		, &	
					v_f				=> this%v_f%v_ptr			, &
					Y				=> this%Y%v_ptr				, &
					Y_int			=> this%Y_int%v_ptr			, &
            		Y_old			=> this%Y_old%v_ptr			, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_droplets	=> this%Y_prod_droplets	, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

			!$omp parallel default(shared)  private(i,j,k,dim,spec,spec_summ,rhs,flux_right,flux_left,lame_coeffs) !, &
			!!$omp& firstprivate(this)
            !!$omp& shared(cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,time_step)
            
			!$omp do collapse(3) schedule(static)						
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
					do spec = 1, species_number

						lame_coeffs		= 1.0_dp				
				  
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
								lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
								lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select						

						Y_rho_int(spec) = 0.5_dp * ( rho_old%cells(i,j,k)*Y_old%pr(spec)%cells(i,j,k) + rho_int%cells(i,j,k)*Y_int%pr(spec)%cells(i,j,k)) 
						
						rhs = 0.0_dp
						do dim = 1, dimensions	
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%CHARM_flux_limiter(rho_int%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
																	 Y_int%pr(spec)%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3)),	&
																	 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
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
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
									flux_left = rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * Y_int%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								else
									flux_left = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
								end if
							end if					

							rhs = rhs - (	flux_right * v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3)&
										-	flux_left * v_f%pr(dim)%cells(dim,i,j,k)* lame_coeffs(dim,1)) / cell_size(1) / lame_coeffs(dim,2)
							
							continue
						end do	
						
						Y_rho_int(spec) = Y_rho_int(spec) + rhs * (0.5_dp * time_step) 
						
						if (this%diffusion_flag)	Y_rho_int(spec) = Y_rho_int(spec) + 0.5_dp * Y_prod_diff%pr(spec)%cells(i,j,k) * time_step
						
						if (this%reactive_flag)		Y_rho_int(spec) = Y_rho_int(spec) + 0.5_dp * Y_prod_chem%pr(spec)%cells(i,j,k) * time_step
						
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
			
			!$omp do collapse(3) schedule(static)
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					spec_summ = 0.0_dp
					do spec = 1,species_number
						spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dp)
					end do
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dp) / spec_summ
					end do
				end if
			end do
			end do
			end do
			!$omp end do

			!$omp end parallel
			
            end associate
			
		!	print *, '5', Y%pr(1)%cells(50,20,1),rho%cells(50,20,1)



    end subroutine	


	subroutine stabilizing_inlet_1D(this,time,stabilized)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time
		logical				,intent(out)	:: stabilized
		
		real(dp)	,dimension(3)		:: cell_size		
		
		real(dp)				:: max_val_lp, left_val, right_val, left_val2, right_val2, left_val3, right_val3, flame_velocity, flame_surface_length, surface_factor
		real(dp)				:: a, b 
		real(dp)				:: time_track, time_delay, time_stabilization
		real(dp), dimension(2), save		:: previous_flame_location = 0.0_dp,	current_flame_location = 0.0_dp

		real(dp), dimension(2), save		:: previous_flame_location_E = 0.0_dp, current_flame_location_E = 0.0_dp, farfield_velocity_E = 0.0_dp
		real(dp), dimension(2), save		:: previous_flame_location_T = 0.0_dp, current_flame_location_T = 0.0_dp, farfield_velocity_T = 0.0_dp
		real(dp), dimension(2), save		:: previous_flame_location_H = 0.0_dp, current_flame_location_H = 0.0_dp, farfield_velocity_H = 0.0_dp

		real(dp), save			:: flame_velocity_E, max_val_E, flame_velocity_T, max_val_T, flame_velocity_H, max_val_H, farfield_velocity = 0.0_dp, av_flame_velocity = 0.0

		real(dp), save			:: previous_time = 0.0_dp, current_time = 0.0_dp, previous_correction_time = 0.0_dp

		real(dp), dimension(20), save	:: flame_velocity_hist = 0.0
		integer		,save			:: track_counter = 0, correction_counter = 0, stabilization_counter = 0, hist_indx

		character(len=200)			:: file_name
		
		integer		,dimension(:), allocatable, save	:: flame_front_index
		real(dp)	,dimension(:), allocatable, save	:: max_val, max_coord
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dp)					:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x, velocity_deviation
		real(dp)	,dimension(2)	:: coords_lp
		integer						:: i_lp, j_lp, lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2
		integer		,save			:: j_lp_E = 1
		
		integer						:: hist_size
		real(dp)	,dimension(20)	:: data_array, test
		real(dp)					:: s, diff, var_s, z
        
		integer :: CO_index, H2O2_index, HO2_index, OH_index, H_index
		integer	:: bound_number,sign
		integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		
		logical	,save			:: correction_flag = .true.
		character(len=20)		:: boundary_type_name
		character(len=20)		:: flame_data_file
        
        integer	,save			:: flame_loc_unit
        integer					:: flame_loc_unit_old, stat
        
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

!		CO_index		= this%chem%chem_ptr%get_chemical_specie_index('CO')
!		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
!		OH_index		= this%chem%chem_ptr%get_chemical_specie_index('OH')
		H_index			= this%chem%chem_ptr%get_chemical_specie_index('H')
!		HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')

		stabilized = .false.
		
		if (.not.allocated(flame_front_index)) then
			allocate(flame_front_index(cons_inner_loop(2,1):cons_inner_loop(2,2)))
			allocate(max_val(cons_inner_loop(2,1):cons_inner_loop(2,2)))
			allocate(max_coord(cons_inner_loop(2,1):cons_inner_loop(2,2)))
        end if
		
        data_array = 0.0_dp
        
		associate (	v				=> this%v%v_ptr				, &
					v_f				=> this%v_f%v_ptr			, &
					T				=> this%T%s_ptr				, &
					Y				=> this%Y%v_ptr				, &
					E_f_prod_chem 		=> this%E_f_prod_chem%s_ptr	, &
					bc				=> this%boundary%bc_ptr)

		time_delay			= 1e-05_dp!1e-05_dp!1e-05_dp!1e-05_dp			
		time_track			= 1e-03_dp!2e-04_dp!1e-05_dp!2e-04_dp
		time_stabilization	= 1e-03_dp!5e-06_dp!1e-05_dp!5e-06_dp		
       
		if ( (time - time_delay) / time_track  > track_counter + 1) then			
					
            if (track_counter == 0) then
 				write(flame_data_file,'(A,I2,A)') 'av_flame_data_',this%load_counter,'.dat'
				open(newunit = flame_loc_unit, file = flame_data_file, status = 'replace', form = 'formatted')

				if (this%load_counter > 1) then
					write(flame_data_file,'(A,I2,A)') 'av_flame_data_',this%load_counter-1,'.dat'
					open(newunit = flame_loc_unit_old, file = flame_data_file, status = 'old', form = 'formatted')
				end if
            end if            
            
			current_time = time
		
			!# 1D front tracer
			flame_front_index	= 0 
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
                
				max_val(j)		= 1.0e07_dp
                max_coord(j)	= cons_inner_loop(1,2) * cell_size(1)
                
				do i = cons_inner_loop(1,1)+10,cons_inner_loop(1,2)-10
					if(bc%bc_markers(i,j,k) == 0) then	
						if ((E_f_prod_chem%cells(i,j,k)) > max_val(j)) then
							max_val(j)				= E_f_prod_chem%cells(i,j,k)
							flame_front_index(j)	= i
						end if
					end if
                end do
                
                if ( flame_front_index(j) /= 0) then
					left_val		= E_f_prod_chem%cells(flame_front_index(j)-1,j,1)
					right_val		= E_f_prod_chem%cells(flame_front_index(j)+1,j,1)	            
              
					max_coord(j)	= (flame_front_index(j) - 0.5_dp)*cell_size(1)
                
					a = (right_val + left_val - 2.0_dp * max_val(j))/2.0_dp/cell_size(1)**2
					b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*max_coord(j) 
                
					max_coord(j) = -b/2.0_dp/a
                end if
			end do
            end do
			
            j_lp		= minloc(max_coord, dim = 1)
            i_lp		= flame_front_index(j_lp)
			max_val_lp	= max_val(j_lp)
			
			j_lp_E		= j_lp

			coords_lp(1)	= (i_lp - 0.5_dp)*cell_size(1)
			coords_lp(2)	= (j_lp - 0.5_dp)*cell_size(1)

			left_val	= E_f_prod_chem%cells(i_lp-1,j_lp,1)
			right_val	= E_f_prod_chem%cells(i_lp+1,j_lp,1)	            

			left_val2	= E_f_prod_chem%cells(i_lp-2,j_lp,1)
			right_val2	= E_f_prod_chem%cells(i_lp+2,j_lp,1)	

			left_val3	= E_f_prod_chem%cells(i_lp-3,j_lp,1)
			right_val3	= E_f_prod_chem%cells(i_lp+3,j_lp,1)            
			a = (right_val + left_val - 2.0_dp * max_val_lp)/2.0_dp/cell_size(1)**2
			b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*coords_lp(1)
			
			max_val_E = max_val_lp

			current_flame_location_E(1) = -b/2.0_dp/a	

			current_flame_location_E(1) = Newton(coords_lp(1),(/left_val3,left_val2,left_val,max_val_E,right_val,right_val2,right_val3/),cell_size(1),1e-10*max_val_E)
			current_flame_location_E(2) = coords_lp(2)
            
			!# 1D front tracer
			flame_front_index	= 0
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
                
				max_val(j) 		= 0.0_dp
                max_coord(j)	= cons_inner_loop(1,2) * cell_size(1)
                
				do i = cons_inner_loop(1,1)+10,cons_inner_loop(1,2)-10
					if(bc%bc_markers(i,j,k) == 0) then	
						if ((T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val(j)) then
							max_val(j) 		= T%cells(i+1,j,k)-T%cells(i-1,j,k)
							flame_front_index(j)	= i
						end if
					end if
                end do                
                
                if ( flame_front_index(j) /= 0) then                
					left_val	= T%cells(i_lp,j_lp,1) - T%cells(i_lp-2,j_lp,1)
					right_val	= T%cells(i_lp+2,j_lp,1) - T%cells(i_lp,j_lp,1)                 
                
					max_coord(j)	= (flame_front_index(j) - 0.5_dp)*cell_size(1)
                
					a = (right_val + left_val - 2.0_dp * max_val(j))/2.0_dp/cell_size(1)**2
					b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*max_coord(j) 
                
					max_coord(j) = -b/2.0_dp/a
                end if
			end do
            end do

            j_lp		= minloc(max_coord, dim = 1)
            i_lp		= flame_front_index(j_lp)
            max_val_lp	= max_val(j_lp)

			coords_lp(1)	= (i_lp - 0.5_dp)*cell_size(1)
			coords_lp(2)	= (j_lp - 0.5_dp)*cell_size(1)

			left_val	= T%cells(i_lp,j_lp,1) - T%cells(i_lp-2,j_lp,1)
			right_val	= T%cells(i_lp+2,j_lp,1) - T%cells(i_lp,j_lp,1)            

			left_val2	= T%cells(i_lp-1,j_lp,1) - T%cells(i_lp-3,j_lp,1)
			right_val2	= T%cells(i_lp+3,j_lp,1) - T%cells(i_lp+1,j_lp,1)

			left_val3	= T%cells(i_lp-2,j_lp,1) - T%cells(i_lp-4,j_lp,1)
			right_val3	= T%cells(i_lp+4,j_lp,1) - T%cells(i_lp+2,j_lp,1)           
            
			a = (right_val + left_val - 2.0_dp * max_val_lp)/2.0_dp/cell_size(1)**2
			b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*coords_lp(1)
			
			max_val_T = max_val_lp

			current_flame_location_T(1) = -b/2.0_dp/a

			current_flame_location_T(1) = Newton(coords_lp(1),(/left_val3,left_val2,left_val,max_val_T,right_val,right_val2,right_val3/),cell_size(1),1e-06*max_val_T)
			current_flame_location_T(2) = coords_lp(2)

			!# 1D front tracer
			flame_front_index	= 0
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
                
				max_val(j) 		= 0.0_dp
                max_coord(j)	= cons_inner_loop(1,2) * cell_size(1)
                
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
					if(bc%bc_markers(i,j,k) == 0) then	
						! max H
						if (abs(Y%pr(H_index)%cells(i,j,k)) > max_val(j)) then
							max_val(j) 		= Y%pr(H_index)%cells(i,j,k)
							flame_front_index(j)	= i
						end if					
					end if
                end do
                
                if ( flame_front_index(j) /= 0) then  
					left_val	= Y%pr(H_index)%cells(i_lp-1,j_lp,1)
					right_val	= Y%pr(H_index)%cells(i_lp+1,j_lp,1)                
                
					max_coord(j)	= (flame_front_index(j) - 0.5_dp)*cell_size(1)
                
					a = (right_val + left_val - 2.0_dp * max_val(j))/2.0_dp/cell_size(1)**2
					b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*max_coord(j) 
                
					max_coord(j) = -b/2.0_dp/a
                end if
			end do
			end do

            j_lp		= minloc(max_coord, dim = 1)
            i_lp		= flame_front_index(j_lp)
			max_val_lp	= max_val(j_lp)
			
			coords_lp(1)	= (i_lp - 0.5_dp)*cell_size(1)
			coords_lp(2)	= (j_lp - 0.5_dp)*cell_size(1)

			left_val	= Y%pr(H_index)%cells(i_lp-1,j_lp,1)
			right_val	= Y%pr(H_index)%cells(i_lp+1,j_lp,1)	            

			left_val2	= Y%pr(H_index)%cells(i_lp-2,j_lp,1)
			right_val2	= Y%pr(H_index)%cells(i_lp+2,j_lp,1)

			left_val3	= Y%pr(H_index)%cells(i_lp-3,j_lp,1)
			right_val3	= Y%pr(H_index)%cells(i_lp+3,j_lp,1)           
            
			a = (right_val + left_val - 2.0_dp * max_val_lp)/2.0_dp/cell_size(1)**2
			b = (right_val - left_val)/2.0_dp/cell_size(1) - 2.0_dp*a*coords_lp(1)
			
			max_val_H = max_val_lp

			current_flame_location_H(1) = -b/2.0_dp/a

			current_flame_location_H(1) = Newton(coords_lp(1),(/left_val3,left_val2,left_val,max_val_H,right_val,right_val2,right_val3/),cell_size(1),1e-06*max_val_H)
			current_flame_location_H(2) = coords_lp(2)

			current_flame_location = current_flame_location_H

			flame_velocity = 0.0_dp
			if(track_counter == 0) then
				farfield_velocity = farfield_velocity_array(1)
				previous_flame_location = current_flame_location
				previous_time = current_time
                if (this%load_counter > 1) then
                    do 
						read (flame_loc_unit_old,'(20E20.12)', iostat = stat) test
                        if (stat == 0) then
                            data_array = test
                        else
                            exit
                        end if
                    end do

					previous_time			= data_array(1)
                    previous_flame_location	= data_array(2:3)
                    av_flame_velocity		= data_array(5)
                    track_counter				= data_array(18)

                    if( track_counter < size(flame_velocity_hist)) then
                        flame_velocity_hist(:mod(track_counter,size(flame_velocity_hist))) = av_flame_velocity
                    else
                        flame_velocity_hist = av_flame_velocity
                    end if
                    
                    flame_velocity_hist(mod(track_counter,size(flame_velocity_hist))+1) = data_array(4)
                end if
 			end if  

        

			if( (track_counter /= 0).and.(track_counter /= nint(data_array(19))).and.(current_flame_location(1) /=  previous_flame_location(1)) )then 
!            if( (correction /= 0).and.(current_flame_location(1) /=  previous_flame_location(1)) )then 
                
				flame_velocity		= (current_flame_location(1) - previous_flame_location(1))/(current_time - previous_time)
                
				flame_velocity_T	= (current_flame_location_T(1) - previous_flame_location_T(1))/(current_time - previous_time)
				flame_velocity_H	= (current_flame_location_H(1) - previous_flame_location_H(1))/(current_time - previous_time)
				flame_velocity_E	= (current_flame_location_E(1) - previous_flame_location_E(1))/(current_time - previous_time)

				flame_velocity_hist(mod(track_counter,size(flame_velocity_hist))+1) = flame_velocity
                
!				av_flame_velocity = sum(flame_velocity_hist)/(count((flame_velocity_hist) > 1e-10_dp) + count((flame_velocity_hist) < -1e-10_dp))         

				av_flame_velocity = sum(flame_velocity_hist)/(count(abs(flame_velocity_hist) > 1e-18_dp)) 
                
                print *, 'Velo hist:'
                print *, flame_velocity_hist
                print *, '----'
                print *, 'Non-zero velocities:', count(abs(flame_velocity_hist) > 1e-18_dp)
                print *, 'Av. flame velocity:', av_flame_velocity
                
                velocity_deviation = 0.0_dp
                do hist_indx = 1, size(flame_velocity_hist)
                    if(abs(flame_velocity_hist(hist_indx)) > 1e-18_dp) then
					!	velocity_deviation = velocity_deviation + (av_flame_velocity - flame_velocity_hist(hist_indx))**2
                        velocity_deviation = velocity_deviation + (av_flame_velocity - flame_velocity_hist(hist_indx))
					end if
                end do
                
                !print *, 'Velocity deviation (before norm):', velocity_deviation

                ! Simplified Mann-Kendall Test (No Tie Adjustment)
                hist_size = size(flame_velocity_hist)
                
                if (count(abs(flame_velocity_hist) > 1e-18_dp) == hist_size) then

					s = 0.0_dp
					! Calculate S statistic (no ties expected)
					do i = 1, hist_size-1
					do j = i+1, hist_size
						diff = flame_velocity_hist(j) - flame_velocity_hist(i)
						if (diff > 0.0_dp) then
							s = s + 1.0_dp    ! Increasing
						else if (diff < 0.0_dp) then
							s = s - 1.0_dp    ! Decreasing
						end if
					end do
					end do

					! Simplified variance (no tie adjustment)
					var_s = real(hist_size)*(hist_size-1)*(2*hist_size+5)/18.0_dp

					! Z-score with continuity correction
					if (s > 0.0_dp) then
						z = (s - 1.0_dp) / sqrt(var_s)
					else if (s < 0.0_dp) then
						z = (s + 1.0_dp) / sqrt(var_s)
					else
						z = 0.0
                    end if
                end if
                
                !velocity_deviation = sqrt(velocity_deviation / (count(abs(flame_velocity_hist) > 1e-18_dp))) / abs(av_flame_velocity) * 100.0_dp
                !velocity_deviation = (velocity_deviation / (count(abs(flame_velocity_hist) > 1e-18_dp))) / abs(av_flame_velocity) * 100.0_dp
                
                !print *, 'Velocity deviation (after norm):', velocity_deviation
                
                !pause

				!if (av_flame_velocity < 0.0_dp) then
				!	correction_flag = .true.
				!else
				!	correction_flag = .false.
				!end if
                
				!if ((track_counter > 10).and.(correction_flag).and.(abs(velocity_deviation) < 1e+01_dp).and.((time - previous_correction_time) / time_stabilization > 1)) then
                if ((count(abs(flame_velocity_hist) > 1e-18_dp) == hist_size).and.(correction_flag))  then
                    print *, abs(z)
                    print *, (time - previous_correction_time) / time_stabilization
                !    pause
					if ((abs(z) < 1.96_dp).and.((time - previous_correction_time) / time_stabilization > 1)) then
						farfield_velocity_array = max(farfield_velocity_array - 0.1_dp*av_flame_velocity, 0.0_dp)
						previous_correction_time = time
						correction_counter = correction_counter + 1 
					end if
                end if
				
                data_array(1)		= time
                data_array(2:3)		= current_flame_location
                data_array(4)		= flame_velocity
                data_array(5)		= av_flame_velocity

                data_array(6)		= farfield_velocity_array(1)
                data_array(7)		= abs(farfield_velocity_array(1))+abs(av_flame_velocity)
                data_array(8)		= max_val_H
                data_array(9:10)	= current_flame_location_T
                data_array(11)		= flame_velocity_T
                data_array(12)		= max_val_T
                data_array(13:14)	= current_flame_location_E
                data_array(15)		= flame_velocity_E
                data_array(16)		= max_val_E
                data_array(17)		= z
                data_array(18)		= correction_counter
                data_array(19)		= stabilization_counter
                data_array(20)		= track_counter
                
                				
				if (mod(track_counter,1) == 0) then
					write (flame_loc_unit,'(20E20.12)') data_array
				end if

				previous_flame_location		= current_flame_location

				previous_flame_location_T	= current_flame_location_T
				previous_flame_location_H	= current_flame_location_H
				previous_flame_location_E	= current_flame_location_E

				previous_time = current_time
				
				if( (track_counter /= 0).and.(abs(av_flame_velocity) < 1e-06))then 
					stabilization_counter = stabilization_counter + 1
				end if
				
			end if

			if(stabilization_counter > 100) then
				stabilized = .true.
			end if
			
			track_counter = track_counter + 1
			
		end if	
			
		end associate

    end subroutine	

	subroutine stabilizing_inlet(this,time)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time
		
		real(dp)	,dimension(3)	:: cell_size		
		
		real(dp)					:: max_grad_temp, left_grad_temp, right_grad_temp, max_CO, left_CO, right_CO, flame_velocity, flame_surface_length, surface_factor
		real(dp)					:: a, b 
		real(dp)					:: time_diff, time_delay
		real(dp), save			:: previous_flame_location = 0.0_dp, current_flame_location = 0.0_dp, farfield_velocity = 0.0_dp
		real(dp), save			:: previous_time = 0.0_dp, current_time = 0.0_dp
		integer		,save			:: correction = 0
		integer						:: flame_front_index
		character(len=200)			:: file_name
		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop

		real(dp)	,dimension(2000,2)	:: lp_coord
		real(dp)	,dimension(2)		:: lp_copy	
		integer(dp)	,dimension(7,2000)	:: chains
		integer			,dimension(7)		:: chain_length		
		
		real(dp)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
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

		time_delay	= 1e-04_dp			
		time_diff	= 5e-05_dp
					
		if ( time > (correction+1)*(time_diff) + time_delay) then			
					
			current_time = time
		
			!# Simple front tracer
			!do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			!do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			!	max_grad_temp = 0.0_dp
			!	max_CO = 0.0_dp
			!	do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
			!		if(bc%bc_markers(i,j,k) == 0) then	
			!		
			!			!! Grad temp
			!			!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_grad_temp) then
			!			!	max_grad_temp = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
			!			!	flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
			!			!	flame_front_index = i
			!			!end if
			!		
			!			! max CO
			!			if (abs(Y%pr(CO_index)%cells(i,j,k)) > max_CO) then
			!				max_CO = Y%pr(CO_index)%cells(i,j,k)
			!				flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
			!				flame_front_index = i
			!			end if
			!	
			!		end if
			!	end do
			!end do
			!end do
			
			!# 2D front tracer
	!		T_flame = 1000.0_dp
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
	!				if( ( abs(lp_coord(lp_index,1) - lp_coord(lp_index2,1)) < 1.0e-08_dp)) then
	!				if( ( abs(lp_coord(lp_index,2) - lp_coord(lp_index2,2)) < 1.0e-08_dp)) then
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
	!		lp_coord(lp_number+1:,:) = -1.0_dp
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
	!		lp_coord(:lp_chain,:)						= -1.0_dp
	!	
	!		lp_chain = chains(maxloc(chain_length,1),maxval(chain_length))
	!		lp_coord(lp_chain:,:)						= -1.0_dp
	!		
	!		
	!!		print *, lp_number, '2'
	!		flame_surface_length = 0.0_dp
	!		current_flame_location = 0.0_dp
	!		lp_number2 = 0
	!		do lp_index = 1, lp_number-1
	!			if((lp_coord(lp_index,1) > 0.0_dp).and.(lp_coord(lp_index+1,1) > 0.0_dp)) then
	!				flame_surface_length = flame_surface_length + sqrt(	(lp_coord(lp_index,1) - lp_coord(lp_index+1,1))**2  + &
	!																	(lp_coord(lp_index,2) - lp_coord(lp_index+1,2))**2)	
	!				current_flame_location = current_flame_location + lp_coord(lp_index,1) 
	!				lp_number2 = lp_number2 + 1
	!			end if
	!		end do
	!		
	!		current_flame_location = current_flame_location/lp_number2
			
			!flame_surface_length = 0.0_dp
			!do j = cons_inner_loop(2,1),cons_inner_loop(2,2)-1 
			!	flame_surface_length = flame_surface_length + sqrt((flame_front_coords(j)-flame_front_coords(j+1))**2 + cell_size(1)**2)
			!end do			
			
			!left_grad_temp	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
			!right_grad_temp	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
   !
			!a = (right_grad_temp + left_grad_temp - 2.0_dp * max_grad_temp)/2.0_dp/cell_size(1)**2
			!b = (max_grad_temp-left_grad_temp)/cell_size(1) - a*(2.0_dp*flame_front_coords(1) - cell_size(1))			
			
			 !left_CO	= Y%pr(CO_index)%cells(flame_front_index-1,1,1)
			 !right_CO	= Y%pr(CO_index)%cells(flame_front_index+1,1,1)		
			
			 !a = (right_CO + left_CO - 2.0_dp * max_CO)/2.0_dp/cell_size(1)**2
			 !b = (max_CO - left_CO)/cell_size(1) - a*(2.0_dp*flame_front_coords(1) - cell_size(1))
			
			 !current_flame_location = -b/2.0_dp/a
			
			!current_flame_location = sum(flame_front_coords) / (cons_inner_loop(2,2) - cons_inner_loop(2,1) + 1)
			!print *, lp_number, '3'
			!current_flame_location = sum(lp_coord(:lp_number,1)) / (lp_number)

			!flame_velocity = 0.0_dp
			!if(correction == 0) then
			!	farfield_velocity = farfield_velocity_array(1)
			!	previous_flame_location = current_flame_location
			!end if

		!	print *, correction, time
		!	pause 
			
			!if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
			!	flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)

		!	print *, flame_velocity, farfield_velocity_array(1), current_flame_location, previous_flame_location,current_time,previous_time
		!	print *, '----'
		!	print *, flame_front_coords
		!	pause
			
			!	farfield_velocity_array = farfield_velocity_array - 0.25_dp*flame_velocity
				
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

			!	write (flame_loc_unit,'(4E14.6)') time, current_flame_location, flame_velocity, farfield_velocity_array(1)! , flame_surface_length, surface_factor
				
			!	previous_flame_location = current_flame_location
			!	previous_time = current_time
				
			!end if

			!correction = correction + 1
			
		end if	
			
		end associate

	end subroutine
	
	subroutine apply_boundary_conditions(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		
		real(dp)	,dimension(3)	:: cell_size		
		
		real(dp)					:: wall_temperature, farfield_temperature, farfield_pressure, farfield_velocity
		
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
											if(v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= T%cells(i,j,k)
											else
												farfield_temperature											= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_temperature()
												T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= farfield_temperature
											end if
										else	
											if(v_f%pr(dim)%cells(dim,i,j,k) < 0.0_dp) then
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
												v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			=  0.0_dp 
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
		real(dp)			,intent(in)	:: time_step
		logical				,intent(in)	:: predictor
		
		real(dp)	:: H_residual
		real(dp)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dp)	,dimension(3)	:: cell_size		
		
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
		real(dp)			,dimension(4)	,intent(in)	:: scalar_array
		real(dp)							,intent(in)	:: velocity
		real(dp)							:: CHARM_flux_limiter
		
		real(dp)	:: B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		
		phi_up	= 0.0_dp
		B_r		= 0.0_dp

		phi_loc = scalar_array(3)	- scalar_array(2)
		if (velocity > 0.0_dp) then
			phi_up = scalar_array(2)	- scalar_array(1)
		else
			phi_up = scalar_array(4)	- scalar_array(3)
		end if
						
		if (phi_loc /= 0) then						
			r	= phi_up / phi_loc
		else
			r	= 0.0_dp
		end if	
  
		!# CHARM
		if ( r > 1.0e-010_dp) then
			s	= 1.0_dp/r
			B_r	= s * (3.0_dp*s + 1.0_dp) / (s+1.0_dp) / (s+1.0_dp)	
		else
			B_r = 0.0_dp
		end if

		!# SUPERBEE
		!B_r	= max(max(0.0_dp,min(2.0_dp*r,1.0_dp)),min(r,2.0_dp))
		
		!# Godunov
		!B_r	= 0.0_dp
        
 		!# Central Difference
		!B_r	= 1.0_dp       

		
		if (velocity > 0.0_dp) then
			CHARM_flux_limiter = scalar_array(2) + B_r * 0.5_dp * phi_up		!# CHARM
		!	CHARM_flux_limiter = scalar_array(2) + B_r * 0.5_dp * phi_loc	!# Others
		else
			CHARM_flux_limiter = scalar_array(3) - B_r * 0.5_dp * phi_up		!# CHARM
		!	CHARM_flux_limiter = scalar_array(3) - B_r * 0.5_dp * phi_loc	!# Others
		end if
		
	end function	

	subroutine calculate_time_step(this)
		class(fds_solver)	,intent(inout)	:: this
		
		real(dp)	:: delta_t_interm, delta_t_interm1, delta_t_interm1x, delta_t_interm2, delta_t_interm3, time_step, time_step2, velocity_value, divergence_value

		real(dp)	,dimension(3)	:: cell_size

		real(dp)	,save	:: time

		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		
		integer	:: sign
		integer :: i,j,k,dim, spec

		time_step		= 10.0_dp !this%initial_time_step
		time_step2		= 1000.0_dp !	
		
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
		
		!$omp parallel default(shared)  private(i,j,k,dim,delta_t_interm,velocity_value) !, &
        !!$omp& firstprivate(this)
		!!$omp& shared(time_step,cons_inner_loop,dimensions,cell_size)

		!$omp do collapse(3) schedule(static) reduction(min:time_step)												

		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dp
				do dim = 1,dimensions
					velocity_value	= velocity_value + abs(v%pr(dim)%cells(i,j,k))/(cell_size(1))       !#L1 norm
				!	velocity_value	= velocity_value + (v%pr(dim)%cells(i,j,k)/(cell_size(1)))**2		!#L2 norm
				end do

				if((velocity_value > 0.0_dp).or.(abs(div_v_int%cells(i,j,k)) > 0.0_dp)) then
					delta_t_interm = 1.0_dp/(velocity_value + abs(div_v_int%cells(i,j,k)))		!# L1 norm
				!	delta_t_interm = 1.0_dp/(sqrt(velocity_value) + abs(div_v_int%cells(i,j,k)))	!# L2 norm

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
							if (D%pr(spec)%cells(i,j,k) > 1e-10_dp) then
								delta_t_interm1x = 1.0_dp/4.0_dp/(D%pr(spec)%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
								if(delta_t_interm1x < delta_t_interm1) delta_t_interm1 = delta_t_interm1x
							end if
						end do
					end if 
					if (this%viscosity_flag)	delta_t_interm2 = 1.0_dp/4.0_dp/(nu%cells(i,j,k)/rho%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
					if (this%heat_trans_flag)	delta_t_interm3 = 1.0_dp/4.0_dp/(kappa%cells(i,j,k)/1000.0_dp/rho%cells(i,j,k))/(dimensions/cell_size(1)/cell_size(1))
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
	!	this%time_step	=	0.25_dp * cell_size(1) / (0.2_dp * sqrt(10*2*Pi))
	!	this%time_step	= 1.0e-04_dp

	!	print *, this%courant_fraction * time_step, time_step2, delta_t_interm1, delta_t_interm2, delta_t_interm3
	
		this%time_step	= min(time_step2,this%courant_fraction * time_step)
		
		this%time_step	= min(5.0e-03_dp, this%time_step)

	!	print *, this%time_step
		
		time = time + this%time_step
		
	!	this%time_step	= time_step
		
		end associate			

	end subroutine

	recursive subroutine V_cycle(this,mesh_iter,time_step,nu_1,nu_2,tol,predictor,V_cycle_depth)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: mesh_iter
		real(dp)			,intent(in)		:: time_step
		integer				,intent(in)		:: nu_1, nu_2
		real(dp)			,intent(in)		:: tol
		logical				,intent(in)		:: predictor
		integer				,intent(in)		:: V_cycle_depth
		
		real(dp), dimension (3,3):: lame_coeffs
		character(len=20)			:: coordinate_system
		
		real(dp)	,dimension(3)	:: cell_size
		integer		,dimension(3,2)	:: cons_utter_loop, cons_inner_loop, subgrid_cons_inner_loop, subgrid_cons_utter_loop, subgrid_cons_inner_loop_higher
		integer		,dimension(3)	:: offset
        real(dp)		:: sub_F_summ, r_summ
		real(dp)		:: beta, farfield_velocity
		
		integer	:: poisson_iteration
		logical	:: converged
		real(dp)	:: a_norm, a_norm_init
		
		character(len=20)		:: boundary_type_name
		
		integer	:: dimensions
		integer	:: factor
		integer	:: plus, sign, bound_number
		integer	:: i,j,k,n,dim,dim1,dim2
		integer	:: ilb, irb, ilt, irt, jlb, jrb, jlt, jrt
        integer	:: ilbr, ilbf, irbr, irbf, iltr, iltf, irtr, irtf
        integer	:: jlbr, jlbf, jrbr, jrbf, jltr, jltf, jrtr, jrtf
        integer	:: lb, rb, lt, rt, ind
        integer :: n_v, n_h, n_d

		integer								    		:: neighbours
        integer,        dimension(:,:),	    allocatable	:: neighbours_indexes
        real(dp),       dimension(:),   	allocatable	:: neighbours_coeffs
        integer,        dimension(:),   	allocatable	:: neighbours_distance
        integer,        dimension(:),   	allocatable	:: neighbours_bound
        integer,        dimension(0:1,0:1,0:1)         	:: neighbours_shifts
        integer,        dimension(:),       allocatable :: nn
        integer,        dimension(3)                    :: shift
        integer										    :: resid
		
		cons_utter_loop	= this%domain%get_local_utter_cells_bounds()	
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		
		subgrid_cons_inner_loop = this%subgrid_cons_inner_loop(mesh_iter,:,:)
		subgrid_cons_utter_loop = this%subgrid_cons_utter_loop(mesh_iter,:,:)
		if (mesh_iter < this%number_of_meshes-1) subgrid_cons_inner_loop_higher = this%subgrid_cons_inner_loop(mesh_iter+1,:,:)
		
		dimensions		= this%domain%get_domain_dimensions()
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		neighbours = 2**dimensions
        allocate(neighbours_indexes(neighbours,3), neighbours_coeffs(neighbours), neighbours_distance(neighbours), neighbours_bound(neighbours))
        allocate(nn(neighbours-1))
        
		factor		= 2**mesh_iter
		
		beta		= 2.0_dp/3.0_dp
		
		cell_size(1)	= cell_size(1)*factor
		
		converged = .false.
		
		coordinate_system	= this%domain%get_coordinate_system_name()

		if (mesh_iter < V_cycle_depth) then	
		
        associate (	ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                    v_f				=> this%v_f%v_ptr			, &
                    v_f_old			=> this%v_f_old%v_ptr		, &            
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

			
            !!$omp& shared(a_norm,beta,factor,nu_1,mesh_iter,coordinate_system,		&
            !!$omp& cons_inner_loop,cons_utter_loop,subgrid_cons_utter_loop,subgrid_cons_inner_loop,subgrid_cons_inner_loop_higher,dimensions,cell_size,converged,predictor,time_step,poisson_iteration,farfield_velocity_array,sub_F_summ,r_summ)	

			poisson_iteration = 0

			do while ((.not.converged).and.(poisson_iteration <= nu_1))
					
				a_norm		= 0.0_dp
				converged	= .false.
                sub_F_summ	= 0.0_dp
                r_summ		= 0.0_dp

                !$omp parallel default(shared)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs) !,			&
			    !!$omp& firstprivate(this)
                
				!$omp do collapse(3) schedule(static) reduction(+:a_norm)
				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then

						lame_coeffs		= 1.0_dp				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select
					
						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dp
								
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
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*sub_R(mesh_iter)%cells(i,j,k)
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k)*lame_coeffs(1,2))

					end if
				end do
				end do
				end do
				!$omp end do

				!$omp do collapse(3) schedule(static)
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
				!if(a_norm/a_norm_init < tolerance) converged = .true.
				
                !$omp end parallel
                
				!if ((poisson_iteration == 0).or.(poisson_iteration == nu_1)) then
				!	print *, mesh_iter, a_norm, poisson_iteration
				!end if				

				poisson_iteration	= poisson_iteration + 1

			end do
			
			!!$omp do collapse(3) schedule(static)
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
			!		lame_coeffs		= 1.0_dp				
		 !
			!		select case(coordinate_system)
			!			case ('cartesian')	
			!				lame_coeffs			= 1.0_dp
			!			case ('cylindrical')
			!				! x -> r, y -> z
			!				lame_coeffs(1,1)	= mesh%mesh(1,offset(1),offset(2),offset(3)) - 0.5_dp*cell_size(1)			
			!				lame_coeffs(1,2)	= mesh%mesh(1,offset(1),offset(2),offset(3))
			!				lame_coeffs(1,3)	= mesh%mesh(1,offset(1),offset(2),offset(3)) + 0.5_dp*cell_size(1)	
			!			case ('spherical')
			!				! x -> r
			!				lame_coeffs(1,1)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)) - 0.5_dp*cell_size(1))**2
			!				lame_coeffs(1,2)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)))**2
			!				lame_coeffs(1,3)	=  (mesh%mesh(1,offset(1),offset(2),offset(3)) + 0.5_dp*cell_size(1))**2
			!		end select
			!	
			!	
			!		sub_R(mesh_iter)%cells(i,j,k) = 0.0_dp
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
			
            !$omp parallel default(shared)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs) !,			&
			!!$omp& firstprivate(this)

			!$omp do collapse(3) schedule(static) reduction(+:sub_F_summ,r_summ)			
			do k = subgrid_cons_inner_loop_higher(3,1),subgrid_cons_inner_loop_higher(3,2)
			do j = subgrid_cons_inner_loop_higher(2,1),subgrid_cons_inner_loop_higher(2,2)
			do i = subgrid_cons_inner_loop_higher(1,1),subgrid_cons_inner_loop_higher(1,2)
				offset = 1
				do dim = 1, dimensions
					offset(dim) = subgrid_cons_inner_loop(dim,1) + 2 * ((i-1) * I_m(dim,1) + (j-1) * I_m(dim,2) + (k-1) * I_m(dim,3))
                end do
				
                if (sub_bc(mesh_iter+1)%cells(i,j,k) == 0) then
                
				!** Piecewise constant restriction (CR)  
					if (dimensions == 1) then
						sub_F(mesh_iter+1)%cells(i,j,k) = 0.5_dp*(sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)))
					end if				

					if (dimensions == 2) then
						sub_F(mesh_iter+1)%cells(i,j,k) = 0.25_dp * (sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)) + sub_R(mesh_iter)%cells(offset(1),offset(2)+1,offset(3)) + sub_R(mesh_iter)%cells(offset(1)+1,offset(2)+1,offset(3)))
                    end if

                    if (dimensions == 3) then
                    	sub_F(mesh_iter+1)%cells(i,j,k) = 0.0625_dp * (	sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3))	+ &
																		sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)) + &
																		sub_R(mesh_iter)%cells(offset(1),offset(2)+1,offset(3)) + &
																		sub_R(mesh_iter)%cells(offset(1),offset(2),offset(3)+1) + &
																		sub_R(mesh_iter)%cells(offset(1)+1,offset(2)+1,offset(3))	+ &
																		sub_R(mesh_iter)%cells(offset(1)+1,offset(2),offset(3)+1)	+ &
																		sub_R(mesh_iter)%cells(offset(1),offset(2)+1,offset(3)+1)	+ &
																		sub_R(mesh_iter)%cells(offset(1)+1,offset(2)+1,offset(3)+1))
					end if

                	select case(coordinate_system)
						case ('cartesian')	
							sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k)
							r_summ = r_summ + 1.0_dp! this%sub_cells_number(mesh_iter+1)
						case ('cylindrical')
							! x -> r, y -> z
							sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k) * sub_mesh(mesh_iter+1)%cells(1,i,j,k)
							r_summ		= r_summ + sub_mesh(mesh_iter+1)%cells(1,i,j,k)
						case ('spherical')
							! x -> r
							sub_F_summ	= sub_F_summ + sub_F(mesh_iter+1)%cells(i,j,k) * sub_mesh(mesh_iter+1)%cells(1,i,j,k)**2
							r_summ		= r_summ + sub_mesh(mesh_iter+1)%cells(1,i,j,k)**2
                    end select
                end if

			end do
			end do
            end do
			!$omp end do
            
            if (this%all_Neumann_flag) then	
				!$omp do collapse(3) schedule(static) 
				do k = subgrid_cons_inner_loop_higher(3,1),subgrid_cons_inner_loop_higher(3,2)
				do j = subgrid_cons_inner_loop_higher(2,1),subgrid_cons_inner_loop_higher(2,2)
				do i = subgrid_cons_inner_loop_higher(1,1),subgrid_cons_inner_loop_higher(1,2)
                    if (sub_bc(mesh_iter+1)%cells(i,j,k) == 0) then
						sub_F(mesh_iter+1)%cells(i,j,k) = sub_F(mesh_iter+1)%cells(i,j,k) - sub_F_summ / r_summ
					end if
				end do
				end do
				end do
				!$omp end do    
            end if
            
			!$omp end parallel				
			
            end associate
                        
			call this%V_cycle(mesh_iter+1,time_step,nu_1,nu_2,tol,predictor,V_cycle_depth)

            associate (	ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                        v_f				=> this%v_f%v_ptr			, &
                        v_f_old			=> this%v_f_old%v_ptr		, &            
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

			!$omp parallel default(shared)  private(i,j,k,n,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs,nn,neighbours_indexes,neighbours_bound,neighbours_coeffs,neighbours_distance,neighbours_shifts,resid,shift) !,			&
			!!$omp& firstprivate(this)
            !!$omp& shared(a_norm,beta,factor,nu_2,mesh_iter,coordinate_system,					&
			!!$omp& cons_inner_loop,cons_utter_loop,subgrid_cons_inner_loop,subgrid_cons_utter_loop,subgrid_cons_inner_loop_higher,dimensions,cell_size,converged,predictor,time_step,farfield_velocity_array,poisson_iteration)

            !$omp do collapse(3) schedule(static)
            do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
			do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
			do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
                    
                sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
                    
                !** translate neighbours into binary number 1 -> (0,0,0); 2 -> (0,0,1); etc
                neighbours_indexes	= 1
                neighbours_coeffs   = 0.0_dp
                neighbours_distance = 0
                shift               = 0
                neighbours_shifts   = 0
                do n = 1, neighbours
                    resid = n - 1
					do dim = dimensions, 1, -1
                            
                        shift(dim) = int(resid / 2**(dim-1)) 
                                                       
						neighbours_indexes(n,dim) = int(resid / 2**(dim-1)) + int(i/2)*I_m(dim,1) + int(j/2)*I_m(dim,2) + int(k/2)*I_m(dim,3)
                        resid = mod(resid, 2**(dim-1))

                        if (( ((i/2) + mod(i,2))*I_m(dim,1) + ((j/2) + mod(j,2))*I_m(dim,2) + ((k/2) + mod(k,2))*I_m(dim,3)) == neighbours_indexes(n,dim)) then
							neighbours_distance(n) = neighbours_distance(n) + 1
                        end if
                    end do
                        
                    neighbours_shifts(shift(1),shift(2),shift(3)) = n
                end do
                    
                    
                neighbours_bound = 0                    
                do n = 1, neighbours
                    bound_number	= sub_bc(mesh_iter+1)%cells(neighbours_indexes(n,1),neighbours_indexes(n,2),neighbours_indexes(n,3))
                            
                    if (( bound_number /= 0 ).and.(neighbours_distance(n) > 1.0e-05)) then
                        neighbours_bound(n) = 0
                        boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
                        select case(boundary_type_name)
                        case('wall')
                            neighbours_bound(n) = neighbours_bound(n) + 1	!# Neumann 
                        case('inlet')
                            neighbours_bound(n) = neighbours_bound(n) + 1	!# Neumann 
                        case('outlet')
                            neighbours_bound(n) = neighbours_bound(n) - 1	!# Dirichlet
                        end select
                    end if
                end do

                do n = 1, neighbours
                    
                    shift = 0
                    if (neighbours_bound(n) == 0) then
                            
                        shift = my_findloc(neighbours_shifts, n) - (/1,1,1/)

                        do dim = 1, dimensions
                            nn(dim)  = neighbours_shifts(   shift(1) * (1 - I_m(dim,1)) + mod(shift(1)+1,2) * I_m(dim,1), &
                                                            shift(2) * (1 - I_m(dim,2)) + mod(shift(2)+1,2) * I_m(dim,2), &
                                                            shift(3) * (1 - I_m(dim,3)) + mod(shift(3)+1,2) * I_m(dim,3))
                        end do
                    
                        !** Piecewise constant prolongation (CP)
                        if (neighbours_distance(n) == dimensions) then
                            neighbours_coeffs(n) = 1.0_dp
                        else
                            neighbours_coeffs(n) = 0.0_dp
                        end if
                        !
                        !** for high-order prolongations there should not be Neumann boundary in the corner cell (cell with highest distance from the fine cell). 
                        
                        if (dimensions == 2) then
                            !** bi linear prolongation (BP) ONLY for 2D
                            !if (neighbours_distance(n) == dimensions) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        !# Neumann boundary
                            !        neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)))*(3 + neighbours_bound(nn(2)))
                            !    else
                            !        !# Dirichlet boundary
                            !        neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)))*(3 + neighbours_bound(nn(2)))
                            !    end if
                            !else if (neighbours_distance(n) == dimensions - 1) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            !# Neumann boundary
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        else
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(1)))
                            !        end if
                            !    else
                            !        !# Dirichlet boundary
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        else
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(1)))
                            !        end if
                            !    end if
                            !else
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        !# Neumann boundary
                            !        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))*(1.0_dp - neighbours_bound(nn(2)))
                            !    else
                            !        !# Dirichlet boundary
                            !        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))*(1.0_dp + neighbours_bound(nn(2)))
                            !    end if
                            !end if
                            !
                            !neighbours_coeffs(n) = 1.0_dp / 16.0_dp * neighbours_coeffs(n)
                        
                        !** Kwak prolongation (KP) ONLY for 2D and 3D
                            !if (neighbours_distance(n) == dimensions) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        !# Neumann boundary
                            !        neighbours_coeffs(n) = (2.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)))
                            !    else
                            !        !# Dirichlet boundary
                            !        neighbours_coeffs(n) = (2.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)))
                            !    end if
                            !else if (neighbours_distance(n) == dimensions - 1) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            !# Neumann boundary
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))
                            !        else
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))
                            !        end if
                            !    else
                            !        !# Dirichlet boundary
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))
                            !        else
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))
                            !        end if
                            !    end if
                            !else
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                            !        !# Neumann boundary
                            !        neighbours_coeffs(n) = 0.0_dp
                            !    else
                            !        !# Dirichlet boundary
                            !        neighbours_coeffs(n) = 0.0_dp
                            !    end if
                            !end if
                            !
                            !
                            !neighbours_coeffs(n) = 1.0_dp / 4.0_dp * neighbours_coeffs(n)
                        end if
                            
                    if (dimensions == 3) then
                            !** tri linear prolongation (BP) ONLY for 3D
                            !if (neighbours_distance(n) == dimensions) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                            !        !# Neumann boundary
                            !        neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !    else
                            !        !# Dirichlet boundary
                            !        neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !    end if
                            !else if (neighbours_distance(n) == dimensions - 1) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                            !        !# Neumann boundary
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !        end if
                            !        if (neighbours_distance(nn(2)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))*(3.0_dp + neighbours_bound(nn(1)))
                            !        end if
                            !        if (neighbours_distance(nn(3)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(3)))*(3.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        end if
                            !    else
                            !        !# Dirichlet boundary
                            !        if (neighbours_distance(nn(1)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !        end if
                            !        if (neighbours_distance(nn(2)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(3)))
                            !        end if
                            !        if (neighbours_distance(nn(3)) == dimensions) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(3)))*(3.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        end if
                            !    end if
                            !else if (neighbours_distance(n) == dimensions - 2) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                            !        !# Neumann boundary
                            !        if ((neighbours_distance(nn(1)) == dimensions - 1).and.(neighbours_distance(nn(2)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))*(1.0_dp - neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !        end if 
                            !        if ((neighbours_distance(nn(3)) == dimensions - 1).and.(neighbours_distance(nn(1)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(3)))*(1.0_dp - neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        end if
                            !        if ((neighbours_distance(nn(2)) == dimensions - 1).and.(neighbours_distance(nn(3)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))*(1.0_dp - neighbours_bound(nn(3)))*(3.0_dp + neighbours_bound(nn(1)))
                            !        end if
                            !    else
                            !        !# Dirichlet boundary
                            !        if ((neighbours_distance(nn(1)) == dimensions - 1).and.(neighbours_distance(nn(2)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))*(1.0_dp + neighbours_bound(nn(2)))*(3.0_dp + neighbours_bound(nn(3)))
                            !        end if 
                            !        if ((neighbours_distance(nn(3)) == dimensions - 1).and.(neighbours_distance(nn(1)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(3)))*(1.0_dp + neighbours_bound(nn(1)))*(3.0_dp + neighbours_bound(nn(2)))
                            !        end if
                            !        if ((neighbours_distance(nn(2)) == dimensions - 1).and.(neighbours_distance(nn(3)) == dimensions - 1)) then
                            !            neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))*(1.0_dp + neighbours_bound(nn(3)))*(3.0_dp + neighbours_bound(nn(1)))
                            !        end if
                            !    end if
                            !else if (neighbours_distance(n) == dimensions - 3) then
                            !    if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                            !    !# Neumann boundary
                            !        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))*(1.0_dp - neighbours_bound(nn(2)))*(1.0_dp - neighbours_bound(nn(3)))
                            !    else
                            !    !# Dirichlet boundary
                            !        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))*(1.0_dp + neighbours_bound(nn(2)))*(1.0_dp + neighbours_bound(nn(3)))
                            !    end if
                            !end if
                            !
                            !neighbours_coeffs(n) = 1.0_dp / 64.0_dp * neighbours_coeffs(n)
                            
                            !** Kwak prolongation (KP) ONLY for 2D and 3D
                            if (neighbours_distance(n) == dimensions) then
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                                    !# Neumann boundary
                                    neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)) + neighbours_bound(nn(3)))
                                else
                                    !# Dirichlet boundary
                                    neighbours_coeffs(n) = (3.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)) + neighbours_bound(nn(3)))
                                end if
                            else if (neighbours_distance(n) == dimensions - 1) then
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                                    !# Neumann boundary
                                    if (neighbours_distance(nn(1)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))
                                    end if
                                    if (neighbours_distance(nn(2)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))
                                    end if
                                    if (neighbours_distance(nn(3)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(3)))
                                    end if
                                else
                                    !# Dirichlet boundary
                                    if (neighbours_distance(nn(1)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))
                                    end if
                                    if (neighbours_distance(nn(2)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))
                                    end if
                                    if (neighbours_distance(nn(3)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(3)))
                                    end if
                                end if
                            else if ((neighbours_distance(n) == dimensions - 3).or.(neighbours_distance(n) == dimensions - 2)) then
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1).or.(neighbours_bound(nn(3)) == 1)) then
                                !# Neumann boundary
                                    neighbours_coeffs(n) = 0.0_dp
                                else
                                !# Dirichlet boundary
                                    neighbours_coeffs(n) = 0.0_dp
                                end if
                            end if
                            
                            neighbours_coeffs(n) = 1.0_dp / 6.0_dp * real(neighbours_coeffs(n), dp)
                            
                        end if
                    end if
                end do
                    
                continue

                do n = 1, neighbours
                    sub_E_old(mesh_iter)%cells(i,j,k) = sub_E_old(mesh_iter)%cells(i,j,k) 	+ neighbours_coeffs(n) * sub_E(mesh_iter+1)%cells(neighbours_indexes(n,1),neighbours_indexes(n,2),neighbours_indexes(n,3))
                end do 
			end do
			end do
			end do				
			!$omp end do
			
            !$omp end parallel
            
            poisson_iteration = 0
		
			do while ((.not.converged).and.(poisson_iteration <= nu_2))
			
				a_norm	= 0.0_dp
				converged = .false.

!$omp parallel default(shared)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs) !,			&
			    !!$omp& firstprivate(this)
                
				!$omp do collapse(3) schedule(static) reduction(+:a_norm)	
				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then
					
						lame_coeffs		= 1.0_dp				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select
					
						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dp
								
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
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*beta*sub_R(mesh_iter)%cells(i,j,k)
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k)*lame_coeffs(1,2))
	 
					end if
				end do
				end do
				end do
				!$omp end do

				!$omp do collapse(3) schedule(static)
				do k = subgrid_cons_utter_loop(3,1),subgrid_cons_utter_loop(3,2)
				do j = subgrid_cons_utter_loop(2,1),subgrid_cons_utter_loop(2,2)
				do i = subgrid_cons_utter_loop(1,1),subgrid_cons_utter_loop(1,2)
					sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
				end do
				end do
				end do	
				!$omp end do
			
				!if(a_norm/a_norm_init < tol) converged = .true.
					
                !$omp end parallel
	
				if ((poisson_iteration == 0).or.(poisson_iteration == nu_2)) then
				!	print *, mesh_iter, a_norm, poisson_iteration
				end if	
				!pause
				!if(a_norm/a_norm_init < tolerance) converged = .true.
					
				poisson_iteration	= poisson_iteration + 1
			end do
			
            end associate
	
		else
			
            associate (	ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
                        v_f				=> this%v_f%v_ptr			, &
                        v_f_old			=> this%v_f_old%v_ptr		, &            
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

			poisson_iteration	= 0
            a_norm_init			= 0.0_dp
		
			do while ((poisson_iteration <= 2000*nu_1).and.(.not.converged))
					
				a_norm		= 0.0_dp
				converged = .false.

       !        !!$omp parallel default(shared)  private(i,j,k,dim,dim1,plus,sign,bound_number,boundary_type_name,farfield_velocity,offset,lame_coeffs,a_norm_init) ,			&
			    !!$omp& firstprivate(this)

                !!$omp do collapse(3) schedule(static) reduction(+:a_norm)			

				do k = subgrid_cons_inner_loop(3,1),subgrid_cons_inner_loop(3,2)
				do j = subgrid_cons_inner_loop(2,1),subgrid_cons_inner_loop(2,2)
				do i = subgrid_cons_inner_loop(1,1),subgrid_cons_inner_loop(1,2)
			
					if(sub_bc(mesh_iter)%cells(i,j,k) == 0) then

						lame_coeffs		= 1.0_dp				
		
						select case(coordinate_system)
							case ('cartesian')	
								lame_coeffs			= 1.0_dp
							case ('cylindrical')
								! x -> r, y -> z
								lame_coeffs(1,1)	= sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1)			
								lame_coeffs(1,2)	= sub_mesh(mesh_iter)%cells(1,i,j,k)
								lame_coeffs(1,3)	= sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1)	
							case ('spherical')
								! x -> r
								lame_coeffs(1,1)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) - 0.5_dp*cell_size(1))**2
								lame_coeffs(1,2)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k))**2
								lame_coeffs(1,3)	=  (sub_mesh(mesh_iter)%cells(1,i,j,k) + 0.5_dp*cell_size(1))**2
						end select

						sub_R(mesh_iter)%cells(i,j,k) = 0.0_dp
								
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
								
						sub_E(mesh_iter)%cells(i,j,k)	= sub_E_old(mesh_iter)%cells(i,j,k) + 1.0_dp/(2.0_dp*dimensions)*sub_R(mesh_iter)%cells(i,j,k) * beta
  
						a_norm = a_norm + abs(sub_R(mesh_iter)%cells(i,j,k)*lame_coeffs(1,2))

					end if
				end do
				end do
				end do
				!!$omp end do
				
				
				!!$omp do collapse(3) schedule(static)	
				do k = subgrid_cons_utter_loop(3,1),subgrid_cons_utter_loop(3,2)
				do j = subgrid_cons_utter_loop(2,1),subgrid_cons_utter_loop(2,2)
				do i = subgrid_cons_utter_loop(1,1),subgrid_cons_utter_loop(1,2)
					sub_E_old(mesh_iter)%cells(i,j,k) = sub_E(mesh_iter)%cells(i,j,k)
				end do
				end do
                end do			
				!!$omp end do
				
                !!$omp end parallel 

				!	print *, a_norm, poisson_iteration
				!	pause
				!print *, sqrt(a_norm), poisson_iteration
				!pause
				!if(a_norm/a_norm_init < tolerance) converged = .true.
				
                if(poisson_iteration == 0) then
                    a_norm_init = a_norm
                end if
                    
				if(a_norm < 1.0e-05) then
                    converged = .true.
                end if
!                print*, a_norm
				poisson_iteration	= poisson_iteration + 1
            end do
		
!            print *, "Deepest level iterations:", poisson_iteration
!            print *, "Deepest level initial and final errors:", a_norm, a_norm_init
		
            continue
		!	print *, mesh_iter, a_norm, poisson_iteration
			
        

        end associate

		end if

		continue
    end subroutine V_cycle
   
    
	subroutine farfield_values_modifier(this, time)
		class(fds_solver)	,intent(in)		:: this
		real(dp)			,intent(in)		:: time
		
		integer				:: j, bound
		character(len=20)	:: boundary_type_name
		real(dp)			:: RND 
		real(dp)	,save	:: farfield_velocity
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
		!		farfield_velocity_array(j) = farfield_velocity + RND*(0.5_dp*farfield_velocity)
		!	end do
		!end if
		
		if (time < 5e-03) then
			farfield_velocity_array = 0.0_dp
		else
			farfield_velocity_array = farfield_velocity
		end if	
	end subroutine
	
	subroutine igniter(this, time)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time
		
		real(dp)				:: duration, delay
		integer					:: i,j,k
		integer	,dimension(3,2)	:: cons_inner_loop
		integer		,save		:: iteration = 0

		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()		
		!delay	= 5.0e-03_dp
		delay	= 0.0e-03_dp
		duration = 800.0e-06_dp
		
		associate (	rho	=> this%rho%s_ptr)
			if ((time <= delay + duration)) then !if ((time > delay).and. (iteration == 0)) then !(time <= delay + duration)) then !(time <= delay + relax)) then
				iteration = iteration + 1
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				!	if ((i-128.5)**2 + (j-40)**2 < 157) then
                    if ((abs(j-256.5) <= 13).and.(i <= 2)) then
						rho%cells(i,j,:)	= 0.1_dp * this%rho_0 !- this%rho_0*18.5_dp/20.0_dp! *(time-delay)/relax
						!T%cells(i,j,:)	= 1500.0_dp
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
		
		real(dp)			:: max_T, min_T, max_grad_T, min_grad_T, lf
		integer				:: max_grad_index, left_min_index, right_min_index
		
		character(len=100)						:: name_string
		character(len=10)						:: fmt
		real(dp)	,dimension(:), allocatable	:: temp

		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		real	,dimension(3)	:: cell_size
		
		integer	:: io_unit
		integer :: table_size
		integer	:: i, spec, dim
		integer	:: ierr
		
		open(newunit = io_unit	, file = trim(table_file), status = 'replace'	, form = 'formatted') 	
		
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

        min_grad_T = 1e06
        do i = cons_inner_loop(1,1),max_grad_index
			if(abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1) < min_grad_T) then
				min_grad_T = abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1)
				left_min_index = i
			end if
        end do

        min_grad_T = 1e06 
        do i = max_grad_index,cons_inner_loop(1,2)
			if(abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1) < min_grad_T) then
				min_grad_T = abs(T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1)
				right_min_index = i
			end if
		end do
        
		min_T = 300.0
		
		lf = (max_T-min_T)/max_grad_T
		
		table_size = max(lf/cell_size(1) - left_min_index, right_min_index - lf/cell_size(1))
 
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
		
		do i = max(cons_inner_loop(1,1), max_grad_index-20*int(lf/cell_size(1))), min(cons_inner_loop(1,2),max_grad_index+20*int(lf/cell_size(1)))
!        do i = max(cons_inner_loop(1,1), max_grad_index-table_size), min(cons_inner_loop(1,2),max_grad_index+table_size)   
!        do i = left_min_index, right_min_index
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
		real(dp)						:: get_time_step
		class(fds_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dp)						:: get_time
		class(fds_solver)	,intent(in)		:: this

		get_time = this%time
    end function
    
    pure function my_findloc(neighbours_shifts, n)
        integer,        dimension(3)                           :: my_findloc
        integer                                 ,intent(in)    :: n
        integer,        dimension(0:1,0:1,0:1)  ,intent(in)    :: neighbours_shifts 
        integer :: i,j,k
        
        do i = 1,2
        do j = 1,2
        do k = 1,2
            if (neighbours_shifts(i-1,j-1,k-1) == n) then
                my_findloc = (/i, j, k/)
            end if
        end do
        end do
        end do
        
    end function
end module
	
!subroutine Weighted_Jacobi_solver_MG(cell_size,max_iter,u,r,f,bc,boundary_types)
!
!	use kind_parameters
!	use global_data
!	use boundary_type_class
!
!    real(dp)								,intent(in)     :: cell_size
!    real(dp)			,dimension(:,:,:)	,intent(inout)  :: u, r, f, bc
!	type(boundary_type)	,dimension(:)		,intent(in)		:: boundary_types
!    integer									,intent(in)     :: max_iter
!
!    real(dp)			,dimension(:,:,:)	,allocatable    :: u_new
!
!	integer	,dimension(3,2)	:: bound_sl
!    integer			:: cells_number(3), dimensions
!    integer			:: i,j,k, dim,iter
!    real(dp)	    :: tolerance, omega, residual_sum, residual_sum_init
!
!	tolerance = 1e-02_dp
!    omega = 2.0_dp/3.0_dp
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
!				r(i,j,k) = r(i,j,k) - (u(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))-2.0_dp*u(i,j,k)+u(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/(cell_size**2)
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
!				u_new(i,j,k) = - (cell_size**2 * f(i,j,k))/(2.0_dp*dimensions)
!			
!				do dim = 1, dimensions
!					u_new(i,j,k) =  u_new(i,j,k) + (u(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + u(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/(2.0_dp*dimensions) 
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
