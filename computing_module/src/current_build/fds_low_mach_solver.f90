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
	use lagrangian_particles_solver_class
	use continuous_particles_solver_class
    use thermal_radiation_solver_class
	
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
    type(timer)     :: fds_radiation_timer
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
		logical			            :: diffusion_flag, viscosity_flag, heat_trans_flag, radiation_flag, reactive_flag, hydrodynamics_flag, CFL_condition_flag, all_Neumann_flag, perturbed_velocity
		real(dp)		            :: courant_fraction
		real(dp)		            :: time, time_step, initial_time_step
        real(dp)    , dimension(3)  :: g
		integer			            :: additional_particles_phases_number
		
		logical			            :: perturbed_velocity_flag, stabilizing_inlet_flag, igniter_flag 
        
		type(viscosity_solver)				:: visc_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(diffusion_solver)				:: diff_solver
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(table_approximated_real_gas)	:: state_eq
        type(thermal_radiation_solver)      :: radiation_solver
        
		type(lagrangian_particles_solver), dimension(:)	    ,allocatable	:: particles_solver			!# Lagrangian particles solver
!		type(continuum_particles_solver), dimension(:)	    ,allocatable	:: particles_solver			!# Continuum particles solver
		
		type(computational_domain)					:: domain
		type(thermophysical_properties_pointer)		:: thermo		
		type(chemical_properties_pointer)			:: chem
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		type(field_scalar_cons_pointer)	:: rho		, rho_int		, rho_old		, T				, T_int			, p				, p_int			, v_s			, mix_mol_mass
		type(field_scalar_cons_pointer)	:: E_f		, E_f_prod_chem	, E_f_prod_heat	, E_f_prod_gd	, E_f_prod_diff	, E_f_prod_rad  , E_f_int		, h_s			, gamma
		type(field_scalar_cons_pointer)	:: p_stat	, p_stat_old	, dp_stat_dt	, p_dyn			, div_v			, div_v_int		, ddiv_v_dt		, H				, H_old			, R
		type(field_scalar_cons_pointer)	:: nu		, kappa
		type(field_scalar_flow_pointer)	:: F_a		, F_b
		
		type(field_vector_cons_pointer)	:: v		, v_prod_gd		, v_prod_visc	, v_prod_sources	, v_int	
		type(field_vector_cons_pointer)	:: Y		, Y_prod_diff	, Y_prod_chem	, Y_int				, Y_old
		type(field_vector_cons_pointer)	:: D
		type(field_vector_flow_pointer)	:: v_f		, v_f_old
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	::  rho_prod_particles, E_f_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  Y_prod_particles, v_prod_particles		
        
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
		procedure	,private	:: eos_corrected_species_face_density
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
	
	type(fds_solver)	function constructor(manager,problem_data_io)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io

		real(dp)						:: calculation_time		
		
		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr
		
		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr		
		
		integer	,dimension(3,2)	:: cons_allocation_bounds		
		integer	,dimension(3,2)	:: cons_utter_loop, cons_inner_loop
		integer	,dimension(3,2)	:: flow_inner_loop, loop
		
		integer				    :: species_number

		type(particles_phase)   :: particles_params
		integer				    :: particles_phase_counter
		
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
		constructor%diffusion_flag		= manager%solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= manager%solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= manager%solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= manager%solver_options%get_chemical_reaction_flag()
        constructor%radiation_flag  	= manager%solver_options%get_thermal_radiation_flag()
		constructor%hydrodynamics_flag	= manager%solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= manager%solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= manager%solver_options%get_CFL_condition_flag()
        
        !# sub solver options
		    constructor%perturbed_velocity_flag	= .false.
            constructor%stabilizing_inlet_flag	= .false.
		    constructor%igniter_flag	        = .false.
            
        constructor%g                   = manager%solver_options%get_grav_acc()
        
		constructor%additional_particles_phases_number	= manager%solver_options%get_additional_particles_phases_number()
		
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
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'specie_mass_fraction')
		constructor%Y%v_ptr						=> vect_ptr%v_ptr
		call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_mass')
		constructor%mix_mol_mass%s_ptr		=> scal_ptr%s_ptr	

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
		call manager%create_vector_field(Y_int		,'specie_mass_fraction_interm'	,'Y_int'	,'chemical')
		constructor%Y_int%v_ptr					=> Y_int
		call manager%create_vector_field(Y_old		,'specie_mass_fraction_old'		,'Y_old'	,'chemical')
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
        
        if (constructor%radiation_flag) then
			constructor%radiation_solver	= thermal_radiation_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_radiation')
			constructor%E_f_prod_rad%s_ptr		=> scal_ptr%s_ptr
		end if

        if(constructor%perturbed_velocity) then
			call manager%create_vector_field(v_prod_sources	,'velocity_production_sources'		,'v_prod_sources'	,'spatial')
			constructor%v_prod_sources%v_ptr	=> v_prod_sources
        end if
        
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
!				constructor%particles_solver(particles_phase_counter)	= continuous_particles_solver_c(manager, particles_params, particles_phase_counter)		!# Continuum particles solver
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_particles(particles_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter						
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)                        !# Lagrangian particles solver
				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr						!# Lagrangian particles solver
    			write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%Y_prod_particles(particles_phase_counter)%v_ptr		=> vect_ptr%v_ptr                
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

			if(constructor%additional_particles_phases_number /= 0) then
				do particles_phase_counter = 1, constructor%additional_particles_phases_number
					call constructor%particles_solver(particles_phase_counter)%set_initial_distributions()
				end do
			end if 			
        end if
        
 		if(constructor%load_counter == 1) then
			call constructor%state_eq%apply_state_equation_for_initial_conditions()
        else
			call constructor%state_eq%apply_state_equation_low_mach_fds(manager%solver_options%get_initial_time_step(),predictor=.true.)
			call constructor%state_eq%apply_boundary_conditions_for_initial_conditions()
        end if	       

        cell_size						= constructor%mesh%mesh_ptr%get_cell_edges_length()
		
		constructor%time				= calculation_time
		constructor%initial_time_step	= manager%solver_options%get_initial_time_step()
		constructor%time_step			= constructor%initial_time_step

		constructor%rho_0				= constructor%rho%s_ptr%cells(cons_inner_loop(1,2),cons_inner_loop(2,2) ,1)!constructor%rho%s_ptr%cells(1,1,1)
!		print *,constructor%rho_0

		species_number = manager%chemistry%chem_ptr%species_number
		
        farfield_velocity = 0.0_dp
		do bound = 1, size(constructor%boundary%bc_ptr%boundary_types)
			boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound)%get_type_name()
			if (boundary_type_name == 'inlet') then
				farfield_velocity = constructor%boundary%bc_ptr%boundary_types(bound)%get_farfield_velocity()
			end if
		end do		

		allocate(farfield_velocity_array(1)) !(cons_allocation_bounds(2,1):cons_allocation_bounds(2,2)))

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
        call manager%create_timer(fds_radiation_timer       ,'FDS radiation solver time'    , 'rad_t')
        call manager%create_timer(fds_viscosity_timer       ,'FDS viscosity solver time'    , 'visc_t')
        call manager%create_timer(fds_multigrid_timer       ,'Multigrid solver time'        , 'mg_t')
 
	end function

	subroutine solve_problem(this,iteration,stop_flag)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: iteration
        logical				,intent(inout)	:: stop_flag
		
		integer	:: particles_phase_counter
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
 
        call fds_radiation_timer%tic()
		if (this%radiation_flag)			call this%radiation_solver%solve_radiation(this%time_step)				
        call fds_radiation_timer%toc(new_iter=.true.)
        
		if (this%perturbed_velocity_flag)		call this%perturb_velocity_field(this%time_step)

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
		end if
         
        
        call fds_timer%toc(new_iter=.true.)
 
		!call this%state_eq%check_conservation_laws()

	end subroutine

	subroutine calculate_interm_Y_predictor(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		
		real(dp)	:: flux_left, flux_right
		real(dp)	:: spec_summ, inv_dx
		
		real(dp)	,dimension(3)	:: cell_size		
		real(dp), dimension (3,3)	:: lame_coeffs
		real(dp), allocatable		:: flux_left_vec(:), flux_right_vec(:), rhs_vec(:), rhoY_new(:)
		
		integer	:: dimensions, species_number, coord_id
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)			:: coordinate_system
		
		integer :: i,j,k,dim,spec,particles_phase_counter
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		inv_dx			= 1.0_dp / cell_size(1)
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		select case(coordinate_system)
			case ('cartesian')
				coord_id = 0
			case ('cylindrical')
				coord_id = 1
			case ('spherical')
				coord_id = 2
			case default
				coord_id = 0
		end select

		associate (	rho				    => this%rho%s_ptr			, &
					rho_int			    => this%rho_int%s_ptr		, &
                    rho_old             => this%rho_old%s_ptr		, &
					v_f				    => this%v_f%v_ptr			, &
					Y				    => this%Y%v_ptr				, &
					Y_int			    => this%Y_int%v_ptr			, &
					Y_prod_diff		    => this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		    => this%Y_prod_chem%v_ptr	, &
					Y_prod_particles	=> this%Y_prod_particles	, &
                    Y_old               => this%Y_old%v_ptr         , &
					mesh			    => this%mesh%mesh_ptr		, &
					bc				    => this%boundary%bc_ptr		, &
					thermo			    => this%thermo%thermo_ptr)
					
			! Optimized scalar predictor:
			!   * species face values are computed as vectors once per face/direction;
			!   * EOS correction is applied once per face, not once per species;
			!   * boundedness and reconstruction are done in the same cell loop;
			!   * coordinate-system string checks are moved out of the inner loop.
			!$omp parallel default(shared) private(i,j,k,dim,spec,spec_summ,flux_right,flux_left,lame_coeffs,flux_left_vec,flux_right_vec,rhs_vec,rhoY_new,particles_phase_counter)
			allocate(flux_left_vec(species_number), flux_right_vec(species_number), rhs_vec(species_number), rhoY_new(species_number))

            !$omp do collapse(3) schedule(static)	
            do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	

					lame_coeffs		= 1.0_dp				
					if (coord_id == 1) then
						! cylindrical: x -> r, y -> z
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
					else if (coord_id == 2) then
						! spherical: x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
					end if

					rhs_vec = 0.0_dp
					do dim = 1, dimensions
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
							call eos_corrected_species_face_vector(rho,Y,thermo%molar_masses,species_number,dim,i,j,k,1, &
								 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),flux_right_vec)
						else
							if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
								do spec = 1,species_number
									flux_right_vec(spec) = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								end do
							else
								do spec = 1,species_number
									flux_right_vec(spec) = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
										Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end do
							end if
						end if
					
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
							call eos_corrected_species_face_vector(rho,Y,thermo%molar_masses,species_number,dim,i,j,k,-1, &
								 v_f%pr(dim)%cells(dim,i,j,k),flux_left_vec)
						else
							if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
								do spec = 1,species_number
									flux_left_vec(spec) = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
										Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end do
							else
								do spec = 1,species_number
									flux_left_vec(spec) = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								end do
							end if
						end if

						do spec = 1,species_number
							rhs_vec(spec) = rhs_vec(spec) -( 	flux_right_vec(spec) * v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) &
														- flux_left_vec(spec)  * v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) * inv_dx / lame_coeffs(dim,2)
						end do
					end do

					do spec = 1, species_number
						rhoY_new(spec) = rho%cells(i,j,k)*Y%pr(spec)%cells(i,j,k) + rhs_vec(spec)*time_step
						if (this%diffusion_flag) rhoY_new(spec) = rhoY_new(spec) + Y_prod_diff%pr(spec)%cells(i,j,k) * time_step
						if (this%reactive_flag)  rhoY_new(spec) = rhoY_new(spec) + Y_prod_chem%pr(spec)%cells(i,j,k) * time_step
						if (this%additional_particles_phases_number /= 0) then
							do particles_phase_counter = 1, this%additional_particles_phases_number
								rhoY_new(spec) = rhoY_new(spec) + Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step
							end do		
						end if						
					end do	

					! Conservative boundedness correction and reconstruction.
					spec_summ = 0.0_dp
					do spec = 1,species_number
						rhoY_new(spec) = max(rhoY_new(spec), 0.0_dp)
						spec_summ = spec_summ + rhoY_new(spec)
					end do
					if (spec_summ > tiny(1.0_dp)) then
						rho_int%cells(i,j,k)	= spec_summ
						do spec = 1, species_number
							Y_int%pr(spec)%cells(i,j,k) = rhoY_new(spec) / rho_int%cells(i,j,k)
						end do
					end if
				end if
			end do
			end do
			end do				
			!$omp end do
			
			! The final copy cannot be merged with the flux-computation loop, because
			! predictor fluxes must see the old rho,Y values in neighbouring cells.
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

			deallocate(flux_left_vec, flux_right_vec, rhs_vec, rhoY_new)
			!$omp end parallel
            end associate
	end subroutine

	subroutine calculate_divergence_v(this,time_step,predictor)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
	
		real(dp)	:: flux_left, flux_right, specie_enthalpy, mixture_cp, dp_dt, D_sum, P_sum, U_sum, div_sum
		
		real(dp)					:: cell_volume
		real(dp)	,dimension(3)	:: cell_size, cell_surface_area
			
		real(dp)					:: energy_source = 1.0e06_dp
		real(dp)	,save			:: time
		integer		,save			:: iter = 1

		real(dp), dimension (3,3)	:: lame_coeffs
		character(len=20)			:: coordinate_system	
		
		integer	:: dimensions, species_number, cells_number
		integer	,dimension(3,2)	:: cons_inner_loop

		integer	:: bound_number, plus, sign
		integer :: i,j,k,dim,spec,particles_phase_counter
		
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
		
		associate (	div_v_int		    => this%div_v_int%s_ptr		, &
					T				    => this%T%s_ptr				, &
					rho				    => this%rho%s_ptr			, &
					p_stat			    => this%p_stat%s_ptr		, &
					gamma			    => this%gamma%s_ptr			, &
					h_s				    => this%h_s%s_ptr			, &
					dp_stat_dt		    => this%dp_stat_dt%s_ptr	, &
					v_f				    => this%v_f%v_ptr			, &
					Y				    => this%Y%v_ptr				, &
					Y_prod_diff		    => this%Y_prod_diff%v_ptr	, &	
					Y_prod_chem		    => this%Y_prod_chem%v_ptr	, &
                    mix_mol_mass        => this%mix_mol_mass%s_ptr  , &
					E_f_prod_chem 	    => this%E_f_prod_chem%s_ptr	, &
                    E_f_prod_rad	    => this%E_f_prod_rad%s_ptr	, &
					E_f_prod_heat	    => this%E_f_prod_heat%s_ptr	, &
					E_f_prod_gd 	    => this%E_f_prod_gd%s_ptr	, &
					E_f_prod_diff	    => this%E_f_prod_diff%s_ptr	, &	
					E_f_prod_particles	=> this%E_f_prod_particles  , &
					thermo			    => this%thermo%thermo_ptr	, &
					mesh			    => this%mesh%mesh_ptr		, &
					bc				    => this%boundary%bc_ptr)

			!$omp parallel default(shared)  private(flux_right,flux_left,i,j,k,dim,spec,mixture_cp,specie_enthalpy,plus,sign,bound_number,cell_volume,cell_surface_area,lame_coeffs) !, &
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
					!if (this%viscosity_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_visc%cells(i,j,k)	![J/m^3/s]
					if (this%heat_trans_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_heat%cells(i,j,k)	![J/m^3/s]
					if (this%diffusion_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_diff%cells(i,j,k)	![J/m^3/s]
					if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_chem%cells(i,j,k)	![J/m^3/s]
                    if (this%radiation_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_rad%cells(i,j,k)	

                    !end if

					!average_molar_mass = 0.0_dp
					!do spec = 1,species_number
					!	average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
					!end do
				 !
					!mix_mol_mass		= 1.0_dp / average_molar_mass

					concs = 0.0_dp
					do spec = 1,species_number
						concs(spec)				= Y%pr(spec)%cells(i,j,k) *  mix_mol_mass%cells(i,j,k) / thermo%molar_masses(spec)
					end do		

					if (this%additional_particles_phases_number /= 0) then
						do particles_phase_counter = 1, this%additional_particles_phases_number
							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + E_f_prod_particles(particles_phase_counter)%s_ptr%cells(i,j,k)	![J/m^3/s]
						end do		
                    end if
              

					do dim = 1,dimensions
					
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
							flux_right = charm_face_value(	rho%cells(i-I_m(dim,1):i+2*I_m(dim,1),j-I_m(dim,2):j+2*I_m(dim,2),k-I_m(dim,3):k+2*I_m(dim,3))*		&
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
							flux_left = charm_face_value(	rho%cells(i-2*I_m(dim,1):i+I_m(dim,1),j-2*I_m(dim,2):j+I_m(dim,2),k-2*I_m(dim,3):k+I_m(dim,3))*		&
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
				
					mixture_cp				= thermo%mixture_cp_molar(T%cells(i,j,k), concs)
					
					div_v_int%cells(i,j,k)	= div_v_int%cells(i,j,k) / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))
					
					do spec = 1, species_number
					
						specie_enthalpy = (thermo%specie_enthalpy_molar(T%cells(i,j,k),spec) - thermo%specie_enthalpy_molar(T_ref,spec))  / thermo%molar_masses(spec)
					
						do dim = 1, dimensions
					
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
								flux_right = this%eos_corrected_species_face_density(rho,Y,dim,i,j,k,1, &
														 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),spec)
							else
								if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
									flux_right = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								else
									flux_right = rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * Y%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end if
							end if	
						
							if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
								flux_left = this%eos_corrected_species_face_density(rho,Y,dim,i,j,k,-1, &
														 v_f%pr(dim)%cells(dim,i,j,k),spec)
							else
								if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
									flux_left = rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								else
									flux_left = rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)
								end if
							end if					

							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mix_mol_mass%cells(i,j,k) / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) *  mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))) * &
																			  (	-  (	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* lame_coeffs(dim,3) * (flux_right	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)	&
																					 -	v_f%pr(dim)%cells(dim,i,j,k)									* lame_coeffs(dim,1) * (flux_left	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)) / lame_coeffs(dim,2))
							continue													 
						end do							
						
						
						if (this%diffusion_flag) 	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mix_mol_mass%cells(i,j,k) / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))) * &
																		  (	Y_prod_diff%pr(spec)%cells(i,j,k))

						if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mix_mol_mass%cells(i,j,k) / thermo%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))) * &
																		  (	Y_prod_chem%pr(spec)%cells(i,j,k))
					end do
					
					D_sum = D_sum + div_v_int%cells(i,j,k) * cell_volume
					P_sum = P_sum + (1.0_dp / p_stat%cells(i,j,k) - 1.0_dp / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k))) * cell_volume 
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
					
					!average_molar_mass = 0.0_dp
					!do spec = 1,species_number
					!	average_molar_mass = average_molar_mass + Y%pr(spec)%cells(i,j,k) / thermo%molar_masses(spec)
     !               end do
     !
					!mix_mol_mass	= 1.0_dp / average_molar_mass				
				
                    concs = 0.0_dp
					do spec = 1,species_number
						concs(spec)				= Y%pr(spec)%cells(i,j,k) *  mix_mol_mass%cells(i,j,k) / thermo%molar_masses(spec)
					end do	    
      
					mixture_cp		= thermo%mixture_cp_molar(T%cells(i,j,k), concs)
				
					dp_stat_dt%cells(i,j,k) = (D_sum - U_sum)/ P_sum

					if (this%all_Neumann_flag) then	
						div_v_int%cells(i,j,k)  = div_v_int%cells(i,j,k) - (1.0_dp / p_stat%cells(i,j,k) - 1.0_dp / (rho%cells(i,j,k) * mixture_cp * T%cells(i,j,k) / mix_mol_mass%cells(i,j,k)))* dp_stat_dt%cells(i,j,k)
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
			real(dp)	:: rho_face
		
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
		
		integer	:: particles_phase_counter
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
            
            associate ( rho_old			    => this%rho_old%s_ptr		    , &
                        rho_int			    => this%rho_int%s_ptr		    , &
                        div_v_int		    => this%div_v_int%s_ptr		    , &
                        ddiv_v_dt		    => this%ddiv_v_dt%s_ptr		    , &
                        F_a				    => this%F_a%s_ptr			    , &
                        vorticity		    => this%vorticity			    , &
                        grad_F_a		    => this%grad_F_a			    , &
                        v_f				    => this%v_f%v_ptr			    , &
                        v_f_old			    => this%v_f_old%v_ptr		    , &
                        v_prod_visc		    => this%v_prod_visc%v_ptr	    , &
                        v_prod_particles    => this%v_prod_particles	    , &
                        v_prod_sources	    => this%v_prod_sources%v_ptr    , &
                        mesh			    => this%mesh%mesh_ptr		    , &
                        bc				    => this%boundary%bc_ptr		)

			!$omp parallel default(shared)  private(i,j,k,dim,dim1,dim2,loop,lame_coeffs,farfield_velocity,rho_face,sign,bound_number,bound_number1,bound_number2,bound_number3,plus,boundary_type_name,boundary_type_name1,boundary_type_name2,boundary_type_name3,bc_coeff) !, &
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
							rho_face = 0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))
							if (abs(this%rho_0 - rho_face) > 1e-010_dp) then
								! Consistent face-centered buoyancy: F contains -((rho-rho0)g)/rho.
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - ((this%rho_0 - rho_face) * this%g(dim)) / rho_face
							end if

							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) * (0.5_dp*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
							if (this%perturbed_velocity)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) * (0.5_dp*(v_prod_sources%pr(dim)%cells(i,j,k) + v_prod_sources%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
						    
                            if (this%additional_particles_phases_number /= 0) then
							    do particles_phase_counter = 1, this%additional_particles_phases_number
								    F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k)  + (1.0_dp/(0.5_dp*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) * (0.5_dp*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))	![m/s^2]						!# Lagrangian particles solver
								    !if (bc%bc_markers(i,j,k) == 0) then
								    !	F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum particles solver
								    !end if
							    end do		
						    end if
                        else
							rho_face = 0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))
							if (abs(this%rho_0 - rho_face) > 1e-010_dp) then
								! Consistent face-centered buoyancy: F contains -((rho-rho0)g)/rho.
								F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) - ((this%rho_0 - rho_face) * this%g(dim)) / rho_face
							end if
								
							if (this%viscosity_flag)		F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) * (0.5_dp*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
							if (this%perturbed_velocity)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) * (0.5_dp*(v_prod_sources%pr(dim)%cells(i,j,k) + v_prod_sources%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
                            if (this%additional_particles_phases_number /= 0) then
							    do particles_phase_counter = 1, this%additional_particles_phases_number
								    F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k)  + (1.0_dp/(0.5_dp*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) * (0.5_dp*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))	![m/s^2]						!# Lagrangian particles solver
								    !if (bc%bc_markers(i,j,k) == 0) then
								    !	F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + 0.5_dp*(v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i,j,k) + v_prod_particles(particles_phase_counter)%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) ![m/s^2] !# Continuum particles solver
								    !end if
							    end do		
						    end if
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

				!$omp parallel default(shared)  private(i,j,k,dim,dim2,loop,plus,sign,bound_number,boundary_type_name,residual,lame_coeffs,farfield_velocity,rho_face) !, &
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
                                            
                                        !    farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                            
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

				tolerance = 1e-03_dp
					
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
                                                    
                                                !    farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                                    
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
												R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		            &
																				- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&
																						+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * 1.0_dp * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
												
                                                !R%cells(i,j,k) = R%cells(i,j,k) + (H_old%cells(i,j,k) - H_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
											
                                            case('inlet')
											!	farfield_velocity = farfield_velocity_array(factor * j)
												farfield_velocity = farfield_velocity_array(1)
                                                
                                            !    farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                                
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
                                                    
                                                !    farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                                    
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
												
                                        !    farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                            
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
	
			!$omp parallel default(shared)  private(i,j,k,dim,dim2,vel_abs,plus,sign,bound_number,boundary_type_name,farfield_density,farfield_pressure) !, &
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

										vel_abs = 0.0_dp
										do dim2 = 1,dimensions
											vel_abs = vel_abs + v%pr(dim2)%cells(i,j,k)**2
										end do
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= &
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*vel_abs  
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

										vel_abs = 0.0_dp
										do dim2 = 1,dimensions
											vel_abs = vel_abs + v%pr(dim2)%cells(i,j,k)**2
										end do
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= &
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*vel_abs   
										if (predictor) then	
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										else
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))*rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										end if                                        
                                        
									case('inlet')
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))

										vel_abs = 0.0_dp
										do dim2 = 1,dimensions
											vel_abs = vel_abs + v%pr(dim2)%cells(i,j,k)**2
										end do
										p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= &
											p_dyn%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) - 0.5_dp*vel_abs   
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
										farfield_velocity		=  farfield_velocity_array(1)
                                        !farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                        
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
		
		real(dp)	:: flux_left, flux_right
		real(dp)	:: spec_summ, inv_dx
		
		real(dp), dimension (3,3):: lame_coeffs
		real(dp)	,dimension(3)	:: cell_size		
		real(dp), allocatable		:: flux_left_vec(:), flux_right_vec(:), rhs_vec(:), rhoY_new(:)
		
		integer	:: dimensions, species_number, coord_id
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)	:: coordinate_system
		
		integer	:: i,j,k,dim,spec,particles_phase_counter

		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		inv_dx			= 1.0_dp / cell_size(1)
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		select case(coordinate_system)
			case ('cartesian')
				coord_id = 0
			case ('cylindrical')
				coord_id = 1
			case ('spherical')
				coord_id = 2
			case default
				coord_id = 0
		end select
		
		associate (	rho				    => this%rho%s_ptr			, &
					rho_int			    => this%rho_int%s_ptr		, &
            		rho_old			    => this%rho_old%s_ptr		, &	
					v_f				    => this%v_f%v_ptr			, &
					Y				    => this%Y%v_ptr				, &
					Y_int			    => this%Y_int%v_ptr			, &
            		Y_old			    => this%Y_old%v_ptr			, &
					Y_prod_diff		    => this%Y_prod_diff%v_ptr	, &
					Y_prod_chem		    => this%Y_prod_chem%v_ptr	, &
					Y_prod_particles	=> this%Y_prod_particles	, &
					mesh			    => this%mesh%mesh_ptr		, &
					bc				    => this%boundary%bc_ptr		, &
					thermo			    => this%thermo%thermo_ptr)

			! Optimized scalar corrector: face-vector EOS correction, integer
			! coordinate selector, and boundedness/reconstruction in the same loop.
			!$omp parallel default(shared) private(i,j,k,dim,spec,spec_summ,flux_right,flux_left,lame_coeffs,flux_left_vec,flux_right_vec,rhs_vec,rhoY_new,particles_phase_counter)
			allocate(flux_left_vec(species_number), flux_right_vec(species_number), rhs_vec(species_number), rhoY_new(species_number))
            
			!$omp do collapse(3) schedule(static)						
			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then	
					lame_coeffs		= 1.0_dp				
					if (coord_id == 1) then
						! cylindrical: x -> r, y -> z
						lame_coeffs(1,1)	=  mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1)			
						lame_coeffs(1,2)	=  mesh%mesh(1,i,j,k)
						lame_coeffs(1,3)	=  mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1)	
					else if (coord_id == 2) then
						! spherical: x -> r
						lame_coeffs(1,1)	=  (mesh%mesh(1,i,j,k) - 0.5_dp*cell_size(1))**2
						lame_coeffs(1,2)	=  (mesh%mesh(1,i,j,k))**2
						lame_coeffs(1,3)	=  (mesh%mesh(1,i,j,k) + 0.5_dp*cell_size(1))**2
					end if

					rhs_vec = 0.0_dp
					do dim = 1, dimensions	
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) < cons_inner_loop(dim,2)) then
							call eos_corrected_species_face_vector(rho_int,Y_int,thermo%molar_masses,species_number,dim,i,j,k,1, &
								 v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)),flux_right_vec)
						else
							if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dp) then
								do spec = 1,species_number
									flux_right_vec(spec) = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
								end do
							else
								do spec = 1,species_number
									flux_right_vec(spec) = rho_int%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * &
										Y_int%pr(spec)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end do
							end if
						end if
					
						if ((i*I_m(dim,1) + j*I_m(dim,2)  + k*I_m(dim,3)) > cons_inner_loop(dim,1)) then
							call eos_corrected_species_face_vector(rho_int,Y_int,thermo%molar_masses,species_number,dim,i,j,k,-1, &
								 v_f%pr(dim)%cells(dim,i,j,k),flux_left_vec)
						else
							if (v_f%pr(dim)%cells(dim,i,j,k) > 0.0_dp) then
								do spec = 1,species_number
									flux_left_vec(spec) = rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) * &
										Y_int%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end do
							else
								do spec = 1,species_number
									flux_left_vec(spec) = rho_int%cells(i,j,k) * Y_int%pr(spec)%cells(i,j,k)
								end do
							end if
						end if					

						do spec = 1,species_number
							rhs_vec(spec) = rhs_vec(spec) - ( 	flux_right_vec(spec) * v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3)&
														- flux_left_vec(spec) * v_f%pr(dim)%cells(dim,i,j,k)* lame_coeffs(dim,1)) * inv_dx / lame_coeffs(dim,2)
						end do
					end do	
					
					do spec = 1, species_number
						rhoY_new(spec) = 0.5_dp * ( rho_old%cells(i,j,k)*Y_old%pr(spec)%cells(i,j,k) + rho_int%cells(i,j,k)*Y_int%pr(spec)%cells(i,j,k)) &
										 + rhs_vec(spec) * (0.5_dp * time_step)
						if (this%diffusion_flag) rhoY_new(spec) = rhoY_new(spec) + 0.5_dp * Y_prod_diff%pr(spec)%cells(i,j,k) * time_step
						if (this%reactive_flag)  rhoY_new(spec) = rhoY_new(spec) + 0.5_dp * Y_prod_chem%pr(spec)%cells(i,j,k) * time_step
						if (this%additional_particles_phases_number /= 0) then
							do particles_phase_counter = 1, this%additional_particles_phases_number
								rhoY_new(spec) = rhoY_new(spec) + Y_prod_particles(particles_phase_counter)%v_ptr%pr(spec)%cells(i,j,k) * time_step
							end do		
						end if									
					end do	

					! Conservative boundedness correction and reconstruction.
					spec_summ = 0.0_dp
					do spec = 1,species_number
						rhoY_new(spec) = max(rhoY_new(spec), 0.0_dp)
						spec_summ = spec_summ + rhoY_new(spec)
					end do
					if (spec_summ > tiny(1.0_dp)) then
						rho%cells(i,j,k) = spec_summ
						do spec = 1,species_number
							Y%pr(spec)%cells(i,j,k) = rhoY_new(spec) / rho%cells(i,j,k)
						end do
					end if
				end if
			end do
			end do
			end do
			!$omp end do

			deallocate(flux_left_vec, flux_right_vec, rhs_vec, rhoY_new)
			!$omp end parallel
			
            end associate
    end subroutine	

	subroutine stabilizing_inlet_1D(this,time,stabilized)
!		! Generalized stabilizing inlet routine with two explicit workflows:
!		!   1) stabilize -> write 1D flamelet -> perturb inlet -> measure S_L -> return stabilized=.true.;
!		!   2) stabilize/anchor flame for long-term observation without flamelet or drift measurement.
!		! Anchored-v3 revision: keeps the heat-release gate and spatial anchor of v2,
!		! but adds a fast capture mode for cases where the initial inlet velocity
!		! convects the flame toward the outlet before the fine controller can react.
!		! The routine keeps the old public signature, but the
!		! controlled observable is no longer a single leading point.  It is the
!		! heat-release-weighted flame centroid projected on the inlet axis.
!		!
!		! Main stabilizing changes relative to the first generalized version:
!		!   1. sample-and-hold control: after every inlet correction the velocity
!		!      history is cleared and the controller waits for inlet ramp/settling;
!		!   2. tracking and control intervals are separated;
!		!   3. the inlet velocity is changed through a target/applied ramp, not by
!		!      an instantaneous boundary jump;
!		!   4. proportional gain is small and adaptively damped after sign changes
!		!      or error growth;
!		!   5. when two measurements bracket V_f=0, the next target is chosen by a
!		!      safeguarded bisection step;
!		!   6. a damped secant step is used only when the measured response slope is
!		!      reliable; otherwise the code falls back to a small proportional step;
!		!   7. the controller now uses a domain-safe-window anchor.  Inside the
!		!      usable domain the control target is V_f only; positional feedback is
!		!      activated only when the flame approaches a guarded boundary.  The
!		!      initial flame coordinate is diagnostic and is not used as an exact
!		!      positional set point;
!		!   8. heat-release validity is checked explicitly.  H-radical and T-gradient
!		!      centroids are diagnostics/fallback locations only and are not allowed
!		!      to drive the feedback after the burning front is lost.
!		!
!		! Convention: front_axis=1 corresponds to the present 1D setup.  For an
!		! inlet normal along another direction, replace this by a boundary-derived
!		! inlet-normal coordinate.
!
		class(fds_solver), intent(inout) :: this
		real(dp),          intent(in)    :: time
		logical,           intent(out)   :: stabilized

		integer, parameter :: max_hist_size = 120
		integer, parameter :: max_diag_hist_size = 200
		integer, parameter :: min_hist_for_control = 12
		integer, parameter :: min_hist_for_capture = 6
		integer, parameter :: stable_required_count = 50
		integer, parameter :: flamelet_required_count = 50
		integer, parameter :: post_flamelet_hold_required_count = 100

		real(dp), parameter :: time_delay_default         = 1.0e-03_dp
		real(dp), parameter :: time_track_default         = 1.0e-04_dp
		real(dp), parameter :: time_control_default       = 5.0e-04_dp
		real(dp), parameter :: time_control_capture       = 2.0e-04_dp
		real(dp), parameter :: response_settle_default    = 5.0e-04_dp
		real(dp), parameter :: response_settle_capture    = 2.0e-04_dp
		real(dp), parameter :: inlet_ramp_time_fine       = 2.0e-03_dp
		real(dp), parameter :: inlet_ramp_time_capture    = 1.0e-03_dp
		real(dp), parameter :: controller_gain_initial    = 2.5e-01_dp
		real(dp), parameter :: controller_gain_min        = 2.0e-02_dp
		real(dp), parameter :: controller_gain_max        = 5.0e-01_dp
		real(dp), parameter :: controller_gain_capture    = 5.0e-01_dp
		real(dp), parameter :: controller_max_fraction    = 5.0e-02_dp
		real(dp), parameter :: controller_max_fraction_capture = 1.0e-01_dp
		real(dp), parameter :: controller_max_fraction_emergency = 2.5e-01_dp
		real(dp), parameter :: min_abs_velocity_step      = 2.0e-05_dp
		real(dp), parameter :: min_abs_velocity_step_capture = 1.0e-03_dp
		real(dp), parameter :: emergency_min_fraction     = 1.0e-01_dp
		real(dp), parameter :: velocity_tolerance_on      = 2.0e-05_dp
		real(dp), parameter :: velocity_tolerance_off     = 5.0e-06_dp
		real(dp), parameter :: filter_alpha               = 1.0e-01_dp
		real(dp), parameter :: secant_relaxation          = 2.0e-01_dp
		real(dp), parameter :: feedback_sign_default      = -1.0_dp
		real(dp), parameter :: heat_release_cut_fraction  = 1.0e-08_dp
		real(dp), parameter :: heat_release_valid_relative = 1.0e-07_dp
		real(dp), parameter :: heat_release_valid_absolute = 1.0e-20_dp
		real(dp), parameter :: position_relaxation_time   = 5.0e-02_dp
		real(dp), parameter :: position_tolerance_cells   = 4.0_dp
		real(dp), parameter :: capture_position_tolerance_cells = 12.0_dp
		real(dp), parameter :: capture_velocity_threshold = 2.0e-02_dp
		real(dp), parameter :: outlet_guard_cells         = 20.0_dp
		real(dp), parameter :: outlet_guard_fraction      = 1.5e-01_dp
		real(dp), parameter :: min_secant_du              = 5.0e-05_dp
		real(dp), parameter :: max_response_slope_abs     = 1.0e+03_dp
		real(dp), parameter :: persistent_error_factor    = 5.0_dp
		real(dp), parameter :: ramp_settle_fraction       = 1.0e-03_dp
		real(dp), parameter :: tiny_weight                = tiny(1.0_dp)
		logical,  parameter :: use_secant_control         = .false.

!		! Workflow selector.
!		! workflow_flamelet_sl:
!		!   Stabilize the flame, write 1D flamelet/chemistry data, perturb the inlet,
!		!   measure the linear drift velocity, write the laminar-burning-velocity result,
!		!   and return stabilized=.true. so the external driver can stop/pause.
!		! workflow_anchor_observation:
!		!   Keep the flame anchored for long-term observation.  No flamelet output, no
!		!   open-loop drift measurement, and no final stop signal is generated.
		integer,  parameter :: workflow_flamelet_sl = 1
		integer,  parameter :: workflow_anchor_observation = 2
		integer,  parameter :: active_workflow = workflow_flamelet_sl
		logical,  parameter :: enable_flamelet_output = (active_workflow == workflow_flamelet_sl)
		logical,  parameter :: enable_drift_measurement = (active_workflow == workflow_flamelet_sl)
		logical,  parameter :: pause_after_sl_measurement = (active_workflow == workflow_flamelet_sl)
		integer,  parameter :: stage_anchor_control = 0
		integer,  parameter :: stage_flamelet_ready = 1
		integer,  parameter :: stage_measurement_ramp = 2
		integer,  parameter :: stage_drift_measurement = 3
		integer,  parameter :: stage_measurement_done = 4
		integer,  parameter :: stage_measurement_failed = 5
		integer,  parameter :: min_hist_for_sl_measurement = 40
		real(dp), parameter :: measurement_delta_fraction = 3.0e-02_dp
		real(dp), parameter :: measurement_delta_sign = 1.0_dp
		real(dp), parameter :: measurement_delta_min = 5.0e-03_dp
		real(dp), parameter :: measurement_delta_max_fraction = 3.0e-01_dp
		real(dp), parameter :: measurement_delta_growth = 2.0_dp
		real(dp), parameter :: measurement_max_duration = 2.0e-02_dp
		real(dp), parameter :: measurement_no_motion_displacement_cells = 0.5_dp
		integer,  parameter :: measurement_max_attempts = 5
		real(dp), parameter :: measurement_ramp_time = 1.0e-03_dp
		real(dp), parameter :: measurement_settle_time = 1.0e-03_dp
		real(dp), parameter :: measurement_min_displacement_cells = 6.0_dp
		real(dp), parameter :: measurement_min_duration = 5.0e-03_dp
		real(dp), parameter :: measurement_r2_min = 9.95e-01_dp
		real(dp), parameter :: measurement_split_slope_rel_tol = 2.5e-01_dp
		real(dp), parameter :: measurement_split_slope_abs_tol = 2.0e-04_dp
		real(dp), parameter :: measurement_residual_cells = 1.0_dp

		integer :: dimensions, species_number, boundary_types
		integer :: front_axis
		integer :: i, j, k, dim, bound_number, specie_number, specie_index
		integer :: H2_index, H_index
		integer :: cons_inner_loop(3,2)
		integer :: active_track_number

		real(dp) :: cell_size(3), cell_volume
		real(dp) :: current_flame_location(3)
		real(dp) :: heat_release_centroid(3)
		real(dp) :: H_centroid(3)
		real(dp) :: Tgrad_centroid(3)
		real(dp) :: current_front_coord
		real(dp) :: flame_velocity_lsq, flame_velocity_filtered
		real(dp) :: diag_flame_velocity_lsq, diag_flame_velocity_filtered
		real(dp) :: control_velocity, position_error, position_velocity, position_control_error
		real(dp) :: front_spread, heat_release_integral
		real(dp) :: heat_release_max, heat_release_valid_limit, H_max, Tgrad_max
		real(dp) :: measured_inlet_velocity, proposed_inlet_velocity
		real(dp) :: du_raw, du_limited, max_velocity_step, X_H2
		real(dp) :: time_delay, time_track, time_control, response_settle_time, inlet_ramp_time
		real(dp) :: time_control_effective, response_settle_effective
		real(dp) :: coord(3), weight, qdot_cut, qdot, hval, tgrad
		real(dp) :: sum_weight, sum_s, sum_s2
		real(dp) :: sum_H_weight, sum_Tgrad_weight
		real(dp) :: s_coord
		real(dp) :: ramp_elapsed, ramp_fraction, ramp_residual, ramp_tolerance
		real(dp) :: target_step_for_log
		real(dp) :: sl_displacement, linear_r2, linear_rms, split_slope_diff
		real(dp) :: measurement_elapsed, measurement_displacement, measurement_delta
		logical :: measurement_linear_ok
		real(dp) :: front_position_tolerance, capture_position_tolerance
		real(dp) :: domain_front_min, domain_front_max, domain_front_length
		real(dp) :: outlet_guard_distance, outlet_distance
		real(dp) :: front_safe_min, front_safe_max
		real(dp), allocatable, save :: time_hist(:), front_coord_hist(:)
		real(dp), allocatable, save :: diag_time_hist(:), diag_front_coord_hist(:)
		real(dp), allocatable, save :: farfield_concentrations(:), concs(:)
		character(len=10), allocatable, save :: farfield_species_names(:)
		character(len=5) :: axis_names(3)
		character(len=20) :: flame_data_file
		character(len=500) :: av_header
		character(len=200), save :: data_table_filename, chem_table_filename
		character(len=100) :: chemical_mechanism
		logical :: found_inlet_farfield, trace_success, flame_detected, control_performed
		logical :: measurement_enabled, ramp_settled, reset_history_after_log
		logical :: capture_mode, emergency_mode

		integer, save :: track_counter = 0
		integer, save :: correction_counter = 0
		integer, save :: stabilization_counter = 0
		integer, save :: hist_count = 0
		integer, save :: diag_hist_count = 0
		integer, save :: same_sign_error_counter = 0
		integer, save :: post_flamelet_hold_counter = 0
		integer, save :: measurement_attempt = 0
		integer, save :: flame_loc_unit = -1
		real(dp), save :: previous_correction_time = -huge(1.0_dp)
		real(dp), save :: filtered_velocity_save = 0.0_dp
		real(dp), save :: diag_filtered_velocity_save = 0.0_dp
		real(dp), save :: previous_front_coord = 0.0_dp
		real(dp), save :: adaptive_gain = 5.0e-02_dp
		real(dp), save :: inlet_velocity_target = 0.0_dp
		real(dp), save :: inlet_velocity_applied = 0.0_dp
		real(dp), save :: ramp_start_time = 0.0_dp
		real(dp), save :: ramp_start_velocity = 0.0_dp
		real(dp), save :: active_inlet_ramp_time = inlet_ramp_time_fine
		real(dp), save :: previous_control_velocity = 0.0_dp
		real(dp), save :: previous_control_inlet_velocity = 0.0_dp
		real(dp), save :: front_reference_coord = 0.0_dp
		real(dp), save :: heat_release_peak_save = 0.0_dp
		real(dp), save :: bracket_u_a = 0.0_dp, bracket_v_a = 0.0_dp
		real(dp), save :: bracket_u_b = 0.0_dp, bracket_v_b = 0.0_dp
		logical, save :: initialized = .false.
		logical, save :: inlet_velocity_initialized = .false.
		logical, save :: output_initialized = .false.
		logical, save :: final_output_written = .false.
		logical, save :: flamelet_output_written = .false.
		logical, save :: sl_output_written = .false.
		logical, save :: have_previous_control_point = .false.
		logical, save :: has_bracket = .false.
		logical, save :: front_reference_initialized = .false.
		integer, save :: control_stage = stage_anchor_control
		real(dp), save :: stabilized_inlet_velocity = 0.0_dp
		real(dp), save :: measurement_inlet_velocity = 0.0_dp
		real(dp), save :: measurement_start_time = 0.0_dp
		real(dp), save :: measurement_start_coord = 0.0_dp
		real(dp), save :: sl_displacement_save = 0.0_dp
		real(dp), save :: measurement_velocity_save = 0.0_dp
		real(dp), save :: measurement_r2_save = 0.0_dp
		real(dp), save :: measurement_rms_save = 0.0_dp
		real(dp), save :: measurement_split_slope_diff_save = 0.0_dp
		real(dp), save :: current_measurement_delta = 0.0_dp

		stabilized = .false.

		dimensions      = this%domain%get_domain_dimensions()
		axis_names      = this%domain%get_axis_names()
		species_number  = this%chem%chem_ptr%species_number
		boundary_types  = this%boundary%bc_ptr%get_boundary_types()
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		front_axis = 1
		if (front_axis > dimensions) front_axis = 1

		cell_volume = 1.0_dp
		do dim = 1, dimensions
			cell_volume = cell_volume * cell_size(dim)
		end do
		front_position_tolerance = max(position_tolerance_cells * cell_size(front_axis), 1.0e-8_dp)
		capture_position_tolerance = max(capture_position_tolerance_cells * cell_size(front_axis), &
			front_position_tolerance)
		domain_front_min = (real(cons_inner_loop(front_axis,1),dp) - 0.5_dp) * cell_size(front_axis)
		domain_front_max = (real(cons_inner_loop(front_axis,2),dp) - 0.5_dp) * cell_size(front_axis)
		domain_front_length = max(domain_front_max - domain_front_min, cell_size(front_axis))
		outlet_guard_distance = max(outlet_guard_cells * cell_size(front_axis), &
			outlet_guard_fraction * domain_front_length)

		H2_index = this%chem%chem_ptr%get_chemical_specie_index('H2')
		H_index  = this%chem%chem_ptr%get_chemical_specie_index('H')

		if (.not. allocated(time_hist)) then
			allocate(time_hist(max_hist_size), front_coord_hist(max_hist_size))
			time_hist        = 0.0_dp
			front_coord_hist = 0.0_dp
		end if

		if (.not. allocated(diag_time_hist)) then
			allocate(diag_time_hist(max_diag_hist_size), diag_front_coord_hist(max_diag_hist_size))
			diag_time_hist        = 0.0_dp
			diag_front_coord_hist = 0.0_dp
		end if

		if (.not. initialized) then
			allocate(concs(species_number))
			concs = 0.0_dp

			found_inlet_farfield = .false.
			do bound_number = 1, boundary_types
				if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
					call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_concentrations(farfield_concentrations)
					call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_species_names(farfield_species_names)
					found_inlet_farfield = .true.
					exit
				end if
			end do

			if (found_inlet_farfield .and. allocated(farfield_species_names)) then
				do specie_number = 1, size(farfield_species_names)
					specie_index = this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
					if (specie_index >= 1 .and. specie_index <= species_number) then
						concs(specie_index) = farfield_concentrations(specie_number)
					end if
				end do
			end if

			if (sum(concs) > 0.0_dp .and. H2_index >= 1 .and. H2_index <= species_number) then
				X_H2 = concs(H2_index) / sum(concs) * 100.0_dp
			else
				X_H2 = 0.0_dp
			end if

			chemical_mechanism = trim(this%chem%chem_ptr%get_chemical_mechanism())
			data_table_filename = 'H2-Air_flamelet_' // trim(chemical_mechanism) // '_' // &
				trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'
			chem_table_filename = 'H2-Air_chem_table_' // trim(chemical_mechanism) // '_' // &
				trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'

			previous_correction_time = time
			adaptive_gain = controller_gain_initial
			initialized = .true.
		end if

		if (.not. inlet_velocity_initialized) then
			inlet_velocity_target  = farfield_velocity_array(1)
			inlet_velocity_applied = farfield_velocity_array(1)
			ramp_start_velocity = inlet_velocity_applied
			ramp_start_time = time
			inlet_velocity_initialized = .true.
		end if

		time_delay           = time_delay_default
		time_track           = time_track_default
		time_control         = time_control_default
		response_settle_time = response_settle_default
		inlet_ramp_time      = active_inlet_ramp_time

!		! Apply a finite-duration linear inlet ramp on every call, not only on
!		! tracking calls.  This avoids an acoustic jump but also avoids the very
!		! long tail of exponential relaxation.
		if (inlet_ramp_time > 0.0_dp) then
			ramp_elapsed = max(time - ramp_start_time, 0.0_dp)
			ramp_fraction = min(ramp_elapsed / inlet_ramp_time, 1.0_dp)
		else
			ramp_fraction = 1.0_dp
		end if
		inlet_velocity_applied = ramp_start_velocity + &
			ramp_fraction * (inlet_velocity_target - ramp_start_velocity)
		farfield_velocity_array(1) = inlet_velocity_applied

		ramp_residual = abs(inlet_velocity_target - inlet_velocity_applied)
		ramp_tolerance = max(ramp_settle_fraction * max(abs(inlet_velocity_target), min_abs_velocity_step), &
			1.0e-12_dp)
		ramp_settled = (ramp_residual <= ramp_tolerance .or. ramp_fraction >= 1.0_dp)

		if ((time - time_delay) / time_track <= real(track_counter + 1, dp)) return

		associate (T => this%T%s_ptr, &
				   Y => this%Y%v_ptr, &
				   E_f_prod_chem => this%E_f_prod_chem%s_ptr, &
				   bc => this%boundary%bc_ptr)

			if (.not. output_initialized) call initialize_output_file()

!			! --- First pass: maxima for robust thresholds and diagnostics.
			heat_release_max = 0.0_dp
			H_max            = 0.0_dp
			Tgrad_max        = 0.0_dp
			do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
				if (bc%bc_markers(i,j,k) /= 0) cycle

				qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
				heat_release_max = max(heat_release_max, qdot)

				if (H_index >= 1 .and. H_index <= species_number) then
					hval = abs(Y%pr(H_index)%cells(i,j,k))
					H_max = max(H_max, hval)
				end if

				tgrad = temperature_gradient_norm(i,j,k)
				Tgrad_max = max(Tgrad_max, tgrad)
			end do
			end do
			end do

			heat_release_peak_save = max(heat_release_peak_save, heat_release_max)
			heat_release_valid_limit = max(heat_release_valid_absolute, &
				heat_release_valid_relative * heat_release_peak_save)

!			! --- Second pass: heat-release centroid and fallback centroids.
			heat_release_centroid = 0.0_dp
			H_centroid            = 0.0_dp
			Tgrad_centroid        = 0.0_dp
			sum_weight            = 0.0_dp
			sum_H_weight          = 0.0_dp
			sum_Tgrad_weight      = 0.0_dp
			sum_s                 = 0.0_dp
			sum_s2                = 0.0_dp
			heat_release_integral = 0.0_dp

			qdot_cut = heat_release_cut_fraction * heat_release_max

			do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
				if (bc%bc_markers(i,j,k) /= 0) cycle

				coord = cell_center_coordinates(i,j,k)
				s_coord = coord(front_axis)

				qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
				if (qdot > qdot_cut) then
					weight = qdot * cell_volume
					heat_release_centroid = heat_release_centroid + weight * coord
					sum_weight = sum_weight + weight
					sum_s  = sum_s  + weight * s_coord
					sum_s2 = sum_s2 + weight * s_coord * s_coord
					heat_release_integral = heat_release_integral + weight
				end if

				if (H_index >= 1 .and. H_index <= species_number) then
					hval = max(Y%pr(H_index)%cells(i,j,k), 0.0_dp)
					if (hval > 0.0_dp) then
						weight = hval * cell_volume
						H_centroid = H_centroid + weight * coord
						sum_H_weight = sum_H_weight + weight
					end if
				end if

				tgrad = temperature_gradient_norm(i,j,k)
				if (tgrad > 0.0_dp) then
					weight = tgrad * cell_volume
					Tgrad_centroid = Tgrad_centroid + weight * coord
					sum_Tgrad_weight = sum_Tgrad_weight + weight
				end if
			end do
			end do
			end do

			trace_success = .false.
			flame_detected = .false.
			current_flame_location = 0.0_dp
			front_spread = 0.0_dp

			if (sum_weight > tiny_weight) then
				current_flame_location = heat_release_centroid / sum_weight
				front_spread = max(sum_s2 / sum_weight - (sum_s / sum_weight)**2, 0.0_dp)
				front_spread = sqrt(front_spread)
				trace_success = .true.
				flame_detected = (heat_release_max >= heat_release_valid_limit)
			else if (sum_H_weight > tiny_weight) then
!				! Diagnostic fallback only.  Do not use this point for feedback control.
				current_flame_location = H_centroid / sum_H_weight
				trace_success = .true.
			else if (sum_Tgrad_weight > tiny_weight) then
!				! Diagnostic fallback only.  Do not use this point for feedback control.
				current_flame_location = Tgrad_centroid / sum_Tgrad_weight
				trace_success = .true.
			end if

			if (.not. trace_success) then
				track_counter = track_counter + 1
				return
			end if

			current_front_coord = current_flame_location(front_axis)
			if (flame_detected .and. .not. front_reference_initialized) then
				front_reference_coord = current_front_coord
				front_reference_initialized = .true.
			end if

			front_safe_min = domain_front_min + outlet_guard_distance
			front_safe_max = domain_front_max - outlet_guard_distance
!			! Do not intersect the safe window with front_reference_coord +/-
!			! capture_position_tolerance.  In coarse-grid runs the flame can settle at a
!			! nearby discrete equilibrium that is slightly shifted from the first detected
!			! centroid.  Treating this harmless offset as a positional error produces
!			! two-level inlet-velocity dithering.  The reference coordinate is therefore
!			! kept for diagnostics, while positional feedback is used only as boundary
!			! protection.
			if (front_safe_min >= front_safe_max) then
				front_safe_min = domain_front_min
				front_safe_max = domain_front_max
			end if

			diag_flame_velocity_lsq = 0.0_dp
			diag_flame_velocity_filtered = 0.0_dp
			if (flame_detected) then
				call append_diagnostic_history(time, current_front_coord)
				diag_flame_velocity_lsq = diagnostic_least_squares_velocity()
				if (diag_hist_count <= 2) then
					diag_flame_velocity_filtered = diag_flame_velocity_lsq
					diag_filtered_velocity_save = diag_flame_velocity_filtered
				else
					diag_flame_velocity_filtered = (1.0_dp - filter_alpha) * diag_filtered_velocity_save + &
						filter_alpha * diag_flame_velocity_lsq
					diag_filtered_velocity_save = diag_flame_velocity_filtered
				end if
			else
				diag_hist_count = 0
				diag_filtered_velocity_save = 0.0_dp
			end if

			control_performed = .false.
			reset_history_after_log = .false.
			du_raw = 0.0_dp
			du_limited = 0.0_dp
			target_step_for_log = 0.0_dp
			flame_velocity_lsq = diag_flame_velocity_lsq
			flame_velocity_filtered = diag_flame_velocity_filtered
			position_error = 0.0_dp
			position_control_error = 0.0_dp
			position_velocity = 0.0_dp
			control_velocity = 0.0_dp
			sl_displacement = sl_displacement_save
			linear_r2 = measurement_r2_save
			linear_rms = measurement_rms_save
			split_slope_diff = measurement_split_slope_diff_save
			measurement_elapsed = 0.0_dp
			measurement_displacement = 0.0_dp
			measurement_linear_ok = .false.
			capture_mode = .false.
			emergency_mode = .false.

			if (front_reference_initialized) then
				position_error = current_front_coord - front_reference_coord
				position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
				position_velocity = position_control_error / max(position_relaxation_time, time_track)
				outlet_distance = domain_front_max - current_front_coord
				capture_mode = (position_control_error /= 0.0_dp) .or. &
					(outlet_distance < outlet_guard_distance)
				emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
					(abs(position_control_error) > capture_position_tolerance)
			end if

			if (capture_mode) then
				response_settle_effective = response_settle_capture
				time_control_effective = time_control_capture
			else
				response_settle_effective = response_settle_time
				time_control_effective = time_control
			end if

!			! Do not collect velocity statistics while the inlet is still ramping,
!			! while the flow is inside the post-correction settling interval, or when
!			! the actual heat-release front is lost.  H/T-gradient fallback points are
!			! logged but never allowed to drive the inlet.
			measurement_enabled = flame_detected .and. front_reference_initialized .and. ramp_settled .and. &
				((control_stage == stage_anchor_control) .or. &
				 (control_stage == stage_flamelet_ready)) .and. &
				((time - previous_correction_time) >= response_settle_effective)

			if (measurement_enabled) then
				call append_front_history(time, current_front_coord)
				flame_velocity_lsq = least_squares_velocity()

				if (hist_count <= 2) then
					flame_velocity_filtered = flame_velocity_lsq
					filtered_velocity_save  = flame_velocity_filtered
				else
					flame_velocity_filtered = (1.0_dp - filter_alpha) * filtered_velocity_save + &
						filter_alpha * flame_velocity_lsq
					filtered_velocity_save = flame_velocity_filtered
				end if

				position_error = current_front_coord - front_reference_coord
				position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
				position_velocity = position_control_error / max(position_relaxation_time, time_track)
				control_velocity = flame_velocity_filtered + position_velocity
				outlet_distance = domain_front_max - current_front_coord
				capture_mode = (position_control_error /= 0.0_dp) .or. &
					(abs(flame_velocity_filtered) > capture_velocity_threshold) .or. &
					(outlet_distance < outlet_guard_distance)
				emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
					(abs(position_control_error) > capture_position_tolerance)
				if (capture_mode) then
					time_control_effective = time_control_capture
				else
					time_control_effective = time_control
				end if
			else if (control_stage == stage_anchor_control) then
				call clear_control_history()
				if (.not. flame_detected) has_bracket = .false.
			else if (.not. flame_detected) then
				call clear_control_history()
				has_bracket = .false.
			end if

!			! --- Open-loop drift measurement after the flame has first been stabilized.
!			!     In this stage the inlet target is fixed.  The flame-front coordinate
!			!     must be well approximated by a straight line in time.  The displacement
!			!     flame-speed estimate is then S_L = U_in - V_f for the current sign
!			!     convention.  For other coordinate conventions only this diagnostic sign
!			!     may need to be changed.
			if (control_stage == stage_measurement_ramp) then
				if (flame_detected .and. ramp_settled .and. &
					(time - ramp_start_time) >= measurement_settle_time) then
					call clear_control_history()
					measurement_start_time = time
					measurement_start_coord = current_front_coord
					control_stage = stage_drift_measurement
				end if
			else if (control_stage == stage_drift_measurement) then
				if (flame_detected .and. ramp_settled) then
					call append_front_history(time, current_front_coord)
					flame_velocity_lsq = least_squares_velocity()
					flame_velocity_filtered = flame_velocity_lsq
					measurement_velocity_save = flame_velocity_lsq
					control_velocity = flame_velocity_lsq
					measurement_elapsed = time - measurement_start_time
					measurement_displacement = current_front_coord - measurement_start_coord
					call drift_linearity_diagnostics(linear_r2, linear_rms, split_slope_diff)
					measurement_r2_save = linear_r2
					measurement_rms_save = linear_rms
					measurement_split_slope_diff_save = split_slope_diff
					sl_displacement = inlet_velocity_target - flame_velocity_lsq
					sl_displacement_save = sl_displacement
					measurement_linear_ok = (hist_count >= min_hist_for_sl_measurement) .and. &
						(measurement_elapsed >= measurement_min_duration) .and. &
						(abs(measurement_displacement) >= measurement_min_displacement_cells * cell_size(front_axis)) .and. &
						(linear_r2 >= measurement_r2_min) .and. &
						(linear_rms <= measurement_residual_cells * cell_size(front_axis)) .and. &
						(split_slope_diff <= max(measurement_split_slope_abs_tol, &
							measurement_split_slope_rel_tol * max(abs(flame_velocity_lsq), velocity_tolerance_on)))
					if (measurement_linear_ok) then
						measurement_velocity_save = flame_velocity_lsq
						control_stage = stage_measurement_done
						stabilization_counter = stable_required_count
						call write_laminar_velocity_once()
					else if (measurement_elapsed >= measurement_max_duration .and. &
						abs(measurement_displacement) < measurement_no_motion_displacement_cells * cell_size(front_axis)) then
!						! The perturbation did not move the heat-release centroid by even a
!						! fraction of a cell.  This usually means numerical pinning or that the
!						! imposed velocity perturbation is too small compared with the discrete
!						! front-location uncertainty.  Retry with a larger open-loop perturbation.
						if (measurement_attempt < measurement_max_attempts) then
							measurement_attempt = measurement_attempt + 1
							current_measurement_delta = min(measurement_delta_growth * max(current_measurement_delta, &
								measurement_delta_min), measurement_delta_max_fraction * &
								max(abs(stabilized_inlet_velocity), min_abs_velocity_step))
							measurement_inlet_velocity = max(stabilized_inlet_velocity + &
								measurement_delta_sign * current_measurement_delta, 0.0_dp)
							inlet_velocity_target = measurement_inlet_velocity
							ramp_start_velocity = inlet_velocity_applied
							ramp_start_time = time
							active_inlet_ramp_time = measurement_ramp_time
							control_stage = stage_measurement_ramp
							sl_displacement_save = 0.0_dp
							measurement_velocity_save = 0.0_dp
							measurement_r2_save = 0.0_dp
							measurement_rms_save = 0.0_dp
							measurement_split_slope_diff_save = 0.0_dp
							call clear_control_history()
						else
							control_stage = stage_measurement_failed
						end if
					end if
				else
					call clear_control_history()
				end if
			end if

!			! --- Feedback control.  The correction is rare and based on the response
!			!     measured after the previous target velocity has settled.
			if (measurement_enabled .and. &
				hist_count >= merge(min_hist_for_capture, min_hist_for_control, capture_mode)) then
				if ((time - previous_correction_time) >= time_control_effective) then
					if (control_action_needed()) then
						measured_inlet_velocity = inlet_velocity_target

						call update_adaptive_gain(control_velocity)
						call update_bracket(measured_inlet_velocity, control_velocity)
						call choose_new_inlet_target(measured_inlet_velocity, control_velocity, &
							proposed_inlet_velocity, du_raw, du_limited)

						if (abs(proposed_inlet_velocity - measured_inlet_velocity) > 0.0_dp) then
							inlet_velocity_target = proposed_inlet_velocity
							target_step_for_log = inlet_velocity_target - measured_inlet_velocity
							if (capture_mode) then
								active_inlet_ramp_time = inlet_ramp_time_capture
							else
								active_inlet_ramp_time = inlet_ramp_time_fine
							end if
							ramp_start_velocity = inlet_velocity_applied
							ramp_start_time = time
							previous_correction_time = time
							correction_counter = correction_counter + 1
							control_performed = .true.
							reset_history_after_log = .true.
						end if

						previous_control_inlet_velocity = measured_inlet_velocity
						previous_control_velocity = control_velocity
						have_previous_control_point = .true.
					end if
				end if
			end if

			if (control_stage == stage_anchor_control) then
				if (measurement_enabled .and. ramp_settled .and. hist_count >= min_hist_for_control .and. &
					abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
					position_control_error == 0.0_dp) then
					stabilization_counter = stabilization_counter + 1
				else if (control_performed .or. .not. flame_detected .or. &
					abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
					position_control_error /= 0.0_dp) then
					stabilization_counter = 0
				end if
			else if (control_stage == stage_flamelet_ready) then
				if (measurement_enabled .and. ramp_settled .and. hist_count >= min_hist_for_control .and. &
					abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
					position_control_error == 0.0_dp) then
					post_flamelet_hold_counter = post_flamelet_hold_counter + 1
					stabilization_counter = post_flamelet_hold_counter
				else if (control_performed .or. .not. flame_detected .or. &
					abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
					position_control_error /= 0.0_dp) then
					post_flamelet_hold_counter = 0
					stabilization_counter = 0
				end if
			else if (control_stage == stage_measurement_done) then
				stabilization_counter = stable_required_count
			end if

			if (control_stage == stage_anchor_control .or. control_stage == stage_flamelet_ready) then
				flame_velocity_lsq = diag_flame_velocity_lsq
				flame_velocity_filtered = diag_flame_velocity_filtered
			else if (control_stage == stage_measurement_done) then
				flame_velocity_lsq = measurement_velocity_save
				flame_velocity_filtered = measurement_velocity_save
				control_velocity = measurement_velocity_save
				sl_displacement = sl_displacement_save
				linear_r2 = measurement_r2_save
				linear_rms = measurement_rms_save
				split_slope_diff = measurement_split_slope_diff_save
			end if

			active_track_number = track_counter
			call write_tracking_line()

			if (reset_history_after_log) call clear_control_history()

			previous_front_coord = current_front_coord
			track_counter = track_counter + 1

			if (control_stage == stage_anchor_control .and. &
				stabilization_counter >= flamelet_required_count) then
				if (active_workflow == workflow_flamelet_sl) then
!					! Scenario 1: the first stabilization stage is complete.  Write the
!					! flamelet before any open-loop velocity perturbation is applied.
					if (enable_flamelet_output) call write_flamelet_tables_once()
					control_stage = stage_flamelet_ready
					stabilization_counter = 0
					post_flamelet_hold_counter = 0
					call clear_control_history()
				else
!					! Scenario 2: keep anchoring indefinitely.  Do not write flamelet
!					! data, do not perturb the inlet, and do not return a stop signal.
					stabilized = .false.
				end if

			else if (control_stage == stage_flamelet_ready .and. &
				post_flamelet_hold_counter >= post_flamelet_hold_required_count) then
				if (enable_drift_measurement .and. dimensions == 1) then
					call start_drift_measurement_ramp(time)
				else if (pause_after_sl_measurement) then
					stabilized = .true.
				end if

			else if (control_stage == stage_measurement_done) then
				if (pause_after_sl_measurement) stabilized = .true.
			end if

		end associate

	contains

		subroutine initialize_output_file()
			integer :: local_dim

			write(flame_data_file,'(A,I0,A)') 'av_flame_data_', this%load_counter, '.dat'
			open(newunit = flame_loc_unit, file = flame_data_file, status = 'replace', form = 'formatted')

			av_header = 'VARIABLES="time" '
			do local_dim = 1, dimensions
				av_header = trim(av_header) // '"xf_' // trim(axis_names(local_dim)) // '" '
			end do
			av_header = trim(av_header) // '"Vfl_lsq" "Vfl_filtered" "Vfl_diag_lsq" "Vfl_diag_filtered" "Vfl_measurement_lsq" "Vcontrol" "pos_error" "x_ref" '
			av_header = trim(av_header) // '"U_in_applied" "U_in_target" '
			av_header = trim(av_header) // '"dU_target" "adaptive_gain" "hist_count" ' // &
				'"measurement_on" "bracket_on" "capture_on" "emergency_on" "flame_detected" '
			av_header = trim(av_header) // '"front_spread" "Qint" "Qmax" "Qvalid" "Hmax" "Tgradmax" '
			av_header = trim(av_header) // '"stage" "SL_disp" "lin_R2" "lin_RMS" "split_dV" '
			av_header = trim(av_header) // '"corr_count" "stab_count" "track_count"'

			write(flame_loc_unit,'(A)') trim(av_header)
			output_initialized = .true.
		end subroutine initialize_output_file

		function cell_center_coordinates(ii,jj,kk) result(xc)
			integer, intent(in) :: ii, jj, kk
			real(dp) :: xc(3)

			xc = 0.0_dp
			xc(1) = (real(ii,dp) - 0.5_dp) * cell_size(1)
			if (dimensions >= 2) xc(2) = (real(jj,dp) - 0.5_dp) * cell_size(2)
			if (dimensions >= 3) xc(3) = (real(kk,dp) - 0.5_dp) * cell_size(3)
		end function cell_center_coordinates

		function temperature_gradient_norm(ii,jj,kk) result(grad_norm)
			integer, intent(in) :: ii, jj, kk
			real(dp) :: grad_norm
			real(dp) :: g2, gd

			g2 = 0.0_dp
			if (dimensions >= 1) then
				if (ii > cons_inner_loop(1,1) .and. ii < cons_inner_loop(1,2)) then
					gd = (this%T%s_ptr%cells(ii+1,jj,kk) - this%T%s_ptr%cells(ii-1,jj,kk)) / (2.0_dp * cell_size(1))
					g2 = g2 + gd * gd
				end if
			end if
			if (dimensions >= 2) then
				if (jj > cons_inner_loop(2,1) .and. jj < cons_inner_loop(2,2)) then
					gd = (this%T%s_ptr%cells(ii,jj+1,kk) - this%T%s_ptr%cells(ii,jj-1,kk)) / (2.0_dp * cell_size(2))
					g2 = g2 + gd * gd
				end if
			end if
			if (dimensions >= 3) then
				if (kk > cons_inner_loop(3,1) .and. kk < cons_inner_loop(3,2)) then
					gd = (this%T%s_ptr%cells(ii,jj,kk+1) - this%T%s_ptr%cells(ii,jj,kk-1)) / (2.0_dp * cell_size(3))
					g2 = g2 + gd * gd
				end if
			end if
			grad_norm = sqrt(g2)
		end function temperature_gradient_norm

		subroutine append_front_history(t_new, s_new)
			real(dp), intent(in) :: t_new, s_new

			if (hist_count < max_hist_size) then
				hist_count = hist_count + 1
				time_hist(hist_count) = t_new
				front_coord_hist(hist_count) = s_new
			else
				time_hist(1:max_hist_size-1) = time_hist(2:max_hist_size)
				front_coord_hist(1:max_hist_size-1) = front_coord_hist(2:max_hist_size)
				time_hist(max_hist_size) = t_new
				front_coord_hist(max_hist_size) = s_new
			end if
		end subroutine append_front_history


		subroutine append_diagnostic_history(t_new, s_new)
			real(dp), intent(in) :: t_new, s_new

			if (diag_hist_count < max_diag_hist_size) then
				diag_hist_count = diag_hist_count + 1
				diag_time_hist(diag_hist_count) = t_new
				diag_front_coord_hist(diag_hist_count) = s_new
			else
				diag_time_hist(1:max_diag_hist_size-1) = diag_time_hist(2:max_diag_hist_size)
				diag_front_coord_hist(1:max_diag_hist_size-1) = diag_front_coord_hist(2:max_diag_hist_size)
				diag_time_hist(max_diag_hist_size) = t_new
				diag_front_coord_hist(max_diag_hist_size) = s_new
			end if
		end subroutine append_diagnostic_history


		subroutine clear_control_history()
			hist_count = 0
			filtered_velocity_save = 0.0_dp
			time_hist = 0.0_dp
			front_coord_hist = 0.0_dp
		end subroutine clear_control_history

		function least_squares_velocity() result(vfit)
			real(dp) :: vfit
			integer :: n
			real(dp) :: t_av, s_av, numerator, denominator

			if (hist_count < 2) then
				vfit = 0.0_dp
				return
			end if

			t_av = sum(time_hist(1:hist_count)) / real(hist_count, dp)
			s_av = sum(front_coord_hist(1:hist_count)) / real(hist_count, dp)
			numerator   = 0.0_dp
			denominator = 0.0_dp
			do n = 1, hist_count
				numerator = numerator + (time_hist(n) - t_av) * (front_coord_hist(n) - s_av)
				denominator = denominator + (time_hist(n) - t_av)**2
			end do

			if (denominator > tiny(denominator)) then
				vfit = numerator / denominator
			else
				vfit = 0.0_dp
			end if
		end function least_squares_velocity



		function diagnostic_least_squares_velocity() result(vfit)
			real(dp) :: vfit
			integer :: n
			real(dp) :: t_av, s_av, numerator, denominator

			if (diag_hist_count < 2) then
				vfit = 0.0_dp
				return
			end if

			t_av = sum(diag_time_hist(1:diag_hist_count)) / real(diag_hist_count, dp)
			s_av = sum(diag_front_coord_hist(1:diag_hist_count)) / real(diag_hist_count, dp)
			numerator   = 0.0_dp
			denominator = 0.0_dp
			do n = 1, diag_hist_count
				numerator = numerator + (diag_time_hist(n) - t_av) * &
					(diag_front_coord_hist(n) - s_av)
				denominator = denominator + (diag_time_hist(n) - t_av)**2
			end do

			if (denominator > tiny(denominator)) then
				vfit = numerator / denominator
			else
				vfit = 0.0_dp
			end if
		end function diagnostic_least_squares_velocity


		subroutine drift_linearity_diagnostics(r2_out, rms_out, split_slope_diff_out)
			real(dp), intent(out) :: r2_out, rms_out, split_slope_diff_out
			integer :: n, mid
			real(dp) :: v_all, t_av, s_av, intercept, ss_tot, ss_res, residual
			real(dp) :: v_first, v_second

			r2_out = 0.0_dp
			rms_out = huge(1.0_dp)
			split_slope_diff_out = huge(1.0_dp)
			if (hist_count < 4) return

			v_all = least_squares_velocity()
			t_av = sum(time_hist(1:hist_count)) / real(hist_count, dp)
			s_av = sum(front_coord_hist(1:hist_count)) / real(hist_count, dp)
			intercept = s_av - v_all * t_av
			ss_tot = 0.0_dp
			ss_res = 0.0_dp
			do n = 1, hist_count
				residual = front_coord_hist(n) - (intercept + v_all * time_hist(n))
				ss_res = ss_res + residual * residual
				ss_tot = ss_tot + (front_coord_hist(n) - s_av)**2
			end do
			rms_out = sqrt(ss_res / real(hist_count, dp))
			if (ss_tot > tiny(ss_tot)) then
				r2_out = max(0.0_dp, 1.0_dp - ss_res / ss_tot)
			else
				r2_out = 0.0_dp
			end if

			mid = hist_count / 2
			v_first = least_squares_velocity_range(1, mid)
			v_second = least_squares_velocity_range(mid + 1, hist_count)
			split_slope_diff_out = abs(v_second - v_first)
		end subroutine drift_linearity_diagnostics

		function least_squares_velocity_range(n_first, n_last) result(vfit)
			integer, intent(in) :: n_first, n_last
			real(dp) :: vfit
			integer :: n, n_local
			real(dp) :: t_av, s_av, numerator, denominator

			n_local = n_last - n_first + 1
			if (n_local < 2) then
				vfit = 0.0_dp
				return
			end if

			t_av = 0.0_dp
			s_av = 0.0_dp
			do n = n_first, n_last
				t_av = t_av + time_hist(n)
				s_av = s_av + front_coord_hist(n)
			end do
			t_av = t_av / real(n_local, dp)
			s_av = s_av / real(n_local, dp)

			numerator = 0.0_dp
			denominator = 0.0_dp
			do n = n_first, n_last
				numerator = numerator + (time_hist(n) - t_av) * (front_coord_hist(n) - s_av)
				denominator = denominator + (time_hist(n) - t_av)**2
			end do
			if (denominator > tiny(denominator)) then
				vfit = numerator / denominator
			else
				vfit = 0.0_dp
			end if
		end function least_squares_velocity_range

		pure function deadband_error(value, tolerance) result(error_out)
			real(dp), intent(in) :: value, tolerance
			real(dp) :: error_out

			if (abs(value) <= tolerance) then
				error_out = 0.0_dp
			else
				error_out = sign(abs(value) - tolerance, value)
			end if
		end function deadband_error

		pure function safe_window_error(value, lower_bound, upper_bound) result(error_out)
			real(dp), intent(in) :: value, lower_bound, upper_bound
			real(dp) :: error_out

			if (value < lower_bound) then
				error_out = value - lower_bound
			else if (value > upper_bound) then
				error_out = value - upper_bound
			else
				error_out = 0.0_dp
			end if
		end function safe_window_error

		logical function control_action_needed()
			control_action_needed = (abs(flame_velocity_filtered) > velocity_tolerance_on) .or. &
				(abs(position_control_error) > 0.0_dp .and. abs(control_velocity) > velocity_tolerance_on)
		end function control_action_needed

		subroutine update_adaptive_gain(v_current)
			real(dp), intent(in) :: v_current

!			The previous version damped the gain whenever the measured velocity
!			magnitude increased.  In the supplied logs this drove the gain to its
!			minimum while the flame was still drifting monotonically.  Here the gain
!			is reduced only after a true sign reversal of the anchored control error.
!			Same-sign errors are treated as persistent drift and the gain is allowed
!			to recover.  The finite step limiter and inlet ramp remain responsible for
!			preventing acoustic forcing.
			if (have_previous_control_point) then
				if (v_current * previous_control_velocity < 0.0_dp) then
					same_sign_error_counter = 0
					adaptive_gain = max(0.7_dp * adaptive_gain, controller_gain_min)
				else
					same_sign_error_counter = same_sign_error_counter + 1
					if (same_sign_error_counter >= 2) then
						adaptive_gain = min(1.10_dp * adaptive_gain, controller_gain_max)
					else
						adaptive_gain = min(1.03_dp * adaptive_gain, controller_gain_max)
					end if
				end if
			else
				same_sign_error_counter = 0
			end if
		end subroutine update_adaptive_gain

		subroutine update_bracket(u_current, v_current)
			real(dp), intent(in) :: u_current, v_current

			if (abs(v_current) <= velocity_tolerance_off) return

			if (have_previous_control_point) then
				if (abs(u_current - previous_control_inlet_velocity) >= min_secant_du .and. &
					v_current * previous_control_velocity < 0.0_dp) then
					bracket_u_a = previous_control_inlet_velocity
					bracket_v_a = previous_control_velocity
					bracket_u_b = u_current
					bracket_v_b = v_current
					has_bracket = .true.
				end if
			end if

			if (has_bracket) then
				if (v_current * bracket_v_a > 0.0_dp) then
					bracket_u_a = u_current
					bracket_v_a = v_current
				else if (v_current * bracket_v_b > 0.0_dp) then
					bracket_u_b = u_current
					bracket_v_b = v_current
				end if
				if (bracket_v_a * bracket_v_b > 0.0_dp) has_bracket = .false.
				if (abs(bracket_u_a - bracket_u_b) < min_secant_du) has_bracket = .false.
			end if
		end subroutine update_bracket

		subroutine choose_new_inlet_target(u_current, v_current, u_new, du_unlimited, du_final)
			real(dp), intent(in)  :: u_current, v_current
			real(dp), intent(out) :: u_new, du_unlimited, du_final
			real(dp) :: dU, dV, response_slope, secant_step, bracket_target
			real(dp) :: gain_effective, max_fraction_effective, min_step_effective

			if (emergency_mode) then
				gain_effective = controller_gain_capture
				max_fraction_effective = controller_max_fraction_emergency
				min_step_effective = max(min_abs_velocity_step_capture, &
					emergency_min_fraction * max(abs(u_current), min_abs_velocity_step))
			else if (capture_mode) then
				gain_effective = controller_gain_capture
				max_fraction_effective = controller_max_fraction_capture
				min_step_effective = min_abs_velocity_step_capture
			else
				gain_effective = adaptive_gain
				max_fraction_effective = controller_max_fraction
				min_step_effective = min_abs_velocity_step
			end if

			max_velocity_step = max(max_fraction_effective * &
				max(abs(u_current), min_abs_velocity_step), min_step_effective)

			du_unlimited = feedback_sign_default * gain_effective * v_current
			if (emergency_mode) then
!				! When the flame is close to leaving the useful domain, the controller
!				! must reduce the inflow even if the filtered velocity estimate is
!				! temporarily small.  The direction still follows feedback_sign_default.
				if (abs(du_unlimited) < min_step_effective) then
					du_unlimited = sign(min_step_effective, feedback_sign_default * max(abs(v_current), velocity_tolerance_on))
				end if
			else if (has_bracket .and. .not. capture_mode) then
				bracket_target = 0.5_dp * (bracket_u_a + bracket_u_b)
				du_unlimited = bracket_target - u_current
			else if (use_secant_control .and. have_previous_control_point .and. .not. capture_mode) then
				dU = u_current - previous_control_inlet_velocity
				dV = v_current - previous_control_velocity
				if (abs(dU) >= min_secant_du .and. abs(dV) > velocity_tolerance_off) then
					response_slope = dV / dU
!					! Accept the secant only if the observed response has the sign required
!					! by the stabilizing feedback.  Otherwise it is most likely contaminated
!					! by delayed flame/acoustic transients and would turn the controller into
!					! positive feedback.
					if (response_slope * feedback_sign_default < 0.0_dp .and. &
						abs(response_slope) > 1.0e-12_dp .and. &
						abs(response_slope) < max_response_slope_abs) then
						secant_step = -secant_relaxation * v_current / response_slope
						if (abs(secant_step) <= 5.0_dp * max_velocity_step) then
							du_unlimited = secant_step
						end if
					end if
				end if
			end if

!			! Even a capture/bisection/secant target is passed through a step limiter.
!			! The applied boundary value is additionally ramped in physical time.
			du_final = min(max(du_unlimited, -max_velocity_step), max_velocity_step)

!			! Persistent nonzero drift must produce a finite control action.  This lower
!			! bound is active only outside the velocity/position deadband.
			if (abs(control_velocity) > persistent_error_factor * velocity_tolerance_on .or. &
				abs(position_control_error) > 0.0_dp) then
				if (abs(du_final) < min_step_effective) then
					if (du_unlimited /= 0.0_dp) then
						du_final = sign(min_step_effective, du_unlimited)
					else
						du_final = feedback_sign_default * sign(min_step_effective, v_current)
					end if
				end if
			end if

			u_new = max(u_current + du_final, 0.0_dp)
		end subroutine choose_new_inlet_target

		subroutine write_flamelet_tables_once()
			if (.not. flamelet_output_written) then
				call this%chem_kin_solver%write_chemical_kinetics_table(chem_table_filename)
				call this%write_data_table(data_table_filename)
				flamelet_output_written = .true.
				final_output_written = .true.
			end if
		end subroutine write_flamelet_tables_once

		subroutine write_laminar_velocity_once()
			integer :: sl_unit
			character(len=200) :: sl_file_name

			if (.not. sl_output_written) then
				write(sl_file_name,'(A,I0,A)') 'laminar_flame_velocity_', this%load_counter, '.dat'
				open(newunit = sl_unit, file = sl_file_name, status = 'replace', form = 'formatted')
				write(sl_unit,'(A)') 'VARIABLES="time" "U_in" "Vfl_measurement" "SL_disp" "lin_R2" "lin_RMS" "split_dV" "attempt"'
				write(sl_unit,'(100E20.12)') time, measurement_inlet_velocity, measurement_velocity_save, &
					sl_displacement_save, measurement_r2_save, measurement_rms_save, &
					measurement_split_slope_diff_save, real(measurement_attempt,dp)
				close(sl_unit)
				sl_output_written = .true.
			end if
		end subroutine write_laminar_velocity_once


		subroutine start_drift_measurement_ramp(t_now)
			real(dp), intent(in) :: t_now

			stabilized_inlet_velocity = inlet_velocity_target
			measurement_attempt = 1
			current_measurement_delta = max(measurement_delta_fraction * &
				max(abs(stabilized_inlet_velocity), min_abs_velocity_step), measurement_delta_min)
			current_measurement_delta = min(current_measurement_delta, measurement_delta_max_fraction * &
				max(abs(stabilized_inlet_velocity), min_abs_velocity_step))
			measurement_inlet_velocity = max(stabilized_inlet_velocity + &
				measurement_delta_sign * current_measurement_delta, 0.0_dp)
			inlet_velocity_target = measurement_inlet_velocity
			ramp_start_velocity = inlet_velocity_applied
			ramp_start_time = t_now
			active_inlet_ramp_time = measurement_ramp_time
			sl_displacement_save = 0.0_dp
			measurement_velocity_save = 0.0_dp
			measurement_r2_save = 0.0_dp
			measurement_rms_save = 0.0_dp
			measurement_split_slope_diff_save = 0.0_dp
			control_stage = stage_measurement_ramp
			stabilization_counter = 0
			post_flamelet_hold_counter = 0
			call clear_control_history()
		end subroutine start_drift_measurement_ramp


		subroutine write_tracking_line()
			real(dp) :: measurement_flag, bracket_flag, capture_flag, emergency_flag, flame_detected_flag

			if (measurement_enabled) then
				measurement_flag = 1.0_dp
			else
				measurement_flag = 0.0_dp
			end if

			if (has_bracket) then
				bracket_flag = 1.0_dp
			else
				bracket_flag = 0.0_dp
			end if

			if (capture_mode) then
				capture_flag = 1.0_dp
			else
				capture_flag = 0.0_dp
			end if

			if (emergency_mode) then
				emergency_flag = 1.0_dp
			else
				emergency_flag = 0.0_dp
			end if

			if (flame_detected) then
				flame_detected_flag = 1.0_dp
			else
				flame_detected_flag = 0.0_dp
			end if

			write(flame_loc_unit,'(100E20.12)') &
				time, current_flame_location(1:dimensions), &
				flame_velocity_lsq, flame_velocity_filtered, diag_flame_velocity_lsq, &
				diag_flame_velocity_filtered, measurement_velocity_save, &
				control_velocity, position_error, front_reference_coord, &
				inlet_velocity_applied, inlet_velocity_target, &
				target_step_for_log, adaptive_gain, real(hist_count,dp), measurement_flag, bracket_flag, &
				capture_flag, emergency_flag, flame_detected_flag, &
				front_spread, heat_release_integral, heat_release_max, heat_release_valid_limit, H_max, Tgrad_max, &
				real(control_stage,dp), sl_displacement, linear_r2, linear_rms, split_slope_diff, &
				real(correction_counter,dp), real(stabilization_counter,dp), real(active_track_number,dp)
		end subroutine write_tracking_line

	end subroutine stabilizing_inlet_1D
	
	
	
	subroutine stabilizing_inlet(this,time)
		class(fds_solver)	,intent(inout)	:: this
		real(dp)			,intent(in)		:: time
		
		real(dp)	,dimension(3)	:: cell_size		
		
		real(dp)					:: max_grad_temp, left_grad_temp, right_grad_temp, max_CO, left_CO, right_CO, flame_velocity, flame_surface_length, surface_factor
		real(dp)					:: a, b 
		real(dp)					:: time_diff, time_delay
		real(dp)	,save			:: previous_flame_location = 0.0_dp, current_flame_location = 0.0_dp, farfield_velocity = 0.0_dp
		real(dp)	,save			:: previous_time = 0.0_dp, current_time = 0.0_dp
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
										if(bc%boundary_types(bound_number)%is_conductive()) then
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
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
								
										h_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			            = h_s%cells(i,j,k) * T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) / T%cells(i,j,k)

										do dim1 = 1, dimensions
											if(dim1 == dim) then
											v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))			    = -v%pr(dim1)%cells(i,j,k)
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
										
										 farfield_velocity		=  farfield_velocity_array(1)
                                        !farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()
                                        
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


	pure function charm_face_value(scalar_array,velocity) result(phi)
		real(dp), dimension(4), intent(in) :: scalar_array
		real(dp), intent(in) :: velocity
		real(dp) :: phi
		real(dp) :: B_r, phi_loc, phi_up, r, s
		
		phi_up = 0.0_dp
		B_r = 0.0_dp
		phi_loc = scalar_array(3) - scalar_array(2)
		if (velocity > 0.0_dp) then
			phi_up = scalar_array(2) - scalar_array(1)
		else
			phi_up = scalar_array(4) - scalar_array(3)
		end if
		if (phi_loc /= 0.0_dp) then
			r = phi_up / phi_loc
		else
			r = 0.0_dp
		end if
		if (r > 1.0e-010_dp) then
			s = 1.0_dp / r
			B_r = s * (3.0_dp*s + 1.0_dp) / (s + 1.0_dp) / (s + 1.0_dp)
		else
			B_r = 0.0_dp
		end if
		if (velocity > 0.0_dp) then
			phi = scalar_array(2) + B_r * 0.5_dp * phi_up
		else
			phi = scalar_array(3) - B_r * 0.5_dp * phi_up
		end if
	end function charm_face_value

	subroutine eos_corrected_species_face_vector(rho_field,Y_field,molar_masses,species_number,dim,i,j,k,face_side,velocity,phi)
		type(field_scalar_cons), intent(in) :: rho_field
		type(field_vector_cons), intent(in) :: Y_field
		real(dp), dimension(:), intent(in) :: molar_masses
		integer, intent(in) :: species_number, dim, i, j, k, face_side
		real(dp), intent(in) :: velocity
		real(dp), dimension(:), intent(out) :: phi

		real(dp), dimension(4) :: scalar_array, rho_over_w_array
		real(dp) :: rhoY_max, rhoY_loc, rho_over_w_face, sum_not_gamma
		integer :: s, n, gamma, offset0, offset, gamma_offset
		integer :: ii, jj, kk

		if (face_side > 0) then
			offset0 = -1
		else
			offset0 = -2
		end if

		if (velocity > 0.0_dp) then
			gamma_offset = offset0 + 1
		else
			gamma_offset = offset0 + 2
		end if

		gamma = 1
		rhoY_max = -huge(1.0_dp)
		ii = i + gamma_offset*I_m(dim,1)
		jj = j + gamma_offset*I_m(dim,2)
		kk = k + gamma_offset*I_m(dim,3)
		do s = 1,species_number
			rhoY_loc = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
			if (rhoY_loc > rhoY_max) then
				rhoY_max = rhoY_loc
				gamma = s
			end if
		end do

		rho_over_w_array = 0.0_dp
		do s = 1,species_number
			do n = 1,4
				offset = offset0 + n - 1
				ii = i + offset*I_m(dim,1)
				jj = j + offset*I_m(dim,2)
				kk = k + offset*I_m(dim,3)
				scalar_array(n) = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
				rho_over_w_array(n) = rho_over_w_array(n) + scalar_array(n) / molar_masses(s)
			end do
			phi(s) = charm_face_value(scalar_array,velocity)
		end do

		! FDS-style EOS preservation: correct the locally dominant carrier species
		! so that sum_k (rho Y_k)_f/W_k = (rho/W)_f at the same CHARM face.
		rho_over_w_face = charm_face_value(rho_over_w_array,velocity)
		sum_not_gamma = 0.0_dp
		do s = 1,species_number
			if (s /= gamma) sum_not_gamma = sum_not_gamma + phi(s) / molar_masses(s)
		end do
		phi(gamma) = molar_masses(gamma) * (rho_over_w_face - sum_not_gamma)
	end subroutine eos_corrected_species_face_vector

	function eos_corrected_species_face_density(this,rho_field,Y_field,dim,i,j,k,face_side,velocity,spec) result(phi)
		class(fds_solver)		,intent(in)	:: this
		type(field_scalar_cons)	,intent(in)	:: rho_field
		type(field_vector_cons)	,intent(in)	:: Y_field
		integer					,intent(in)	:: dim, i, j, k, face_side, spec
		real(dp)				,intent(in)	:: velocity
		real(dp)						:: phi

		! Fast EOS-preserving face interpolation.
		!
		! The previous implementation allocated a temporary rhoY_face(:) array and
		! recomputed the CHARM-limited face value for every species every time this
		! function was called. Because the caller invokes the routine inside the
		! species loop, this produced O(Nspec^2) work plus allocation/deallocation
		! overhead for every face. Here we keep the same correction idea, but avoid
		! dynamic allocation and compute the full multi-species correction only for
		! the locally dominant carrier species. For all other species, the routine
		! returns the ordinary CHARM-limited conservative density rho*Y_k.
		!
		! The corrected species is chosen from the upwind-adjacent cell at the face.
		! In hydrogen-air kernels it is practically always N2, so this is equivalent
		! to the FDS most-abundant-species correction while being much cheaper.
		real(dp), dimension(4)	:: scalar_array, rho_over_w_array
		real(dp)				:: rho_over_w_face, sum_not_gamma, rhoY_max, rhoY_loc
		integer					:: species_number, n, s, gamma, offset0, offset, gamma_offset
		integer					:: ii, jj, kk

		species_number = this%chem%chem_ptr%species_number

		if (face_side > 0) then
			offset0 = -1
		else
			offset0 = -2
		end if

		! For the four-point stencil [upstream-upstream, upstream, downstream,
		! downstream-downstream], CHARM uses element 2 for positive velocity and
		! element 3 for negative velocity. Use that same cell to pick the corrected
		! carrier species.
		if (velocity > 0.0_dp) then
			gamma_offset = offset0 + 1
		else
			gamma_offset = offset0 + 2
		end if

		gamma = 1
		rhoY_max = -huge(1.0_dp)
		ii = i + gamma_offset*I_m(dim,1)
		jj = j + gamma_offset*I_m(dim,2)
		kk = k + gamma_offset*I_m(dim,3)
		do s = 1,species_number
			rhoY_loc = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
			if (rhoY_loc > rhoY_max) then
				rhoY_max = rhoY_loc
				gamma = s
			end if
		end do

		if (spec /= gamma) then
			do n = 1,4
				offset = offset0 + n - 1
				ii = i + offset*I_m(dim,1)
				jj = j + offset*I_m(dim,2)
				kk = k + offset*I_m(dim,3)
				scalar_array(n) = rho_field%cells(ii,jj,kk) * Y_field%pr(spec)%cells(ii,jj,kk)
			end do
			phi = charm_face_value(scalar_array,velocity)
			return
		end if

		! Only the dominant carrier species absorbs the EOS correction.
		rho_over_w_array = 0.0_dp
		sum_not_gamma = 0.0_dp
		do s = 1,species_number
			do n = 1,4
				offset = offset0 + n - 1
				ii = i + offset*I_m(dim,1)
				jj = j + offset*I_m(dim,2)
				kk = k + offset*I_m(dim,3)
				scalar_array(n) = rho_field%cells(ii,jj,kk) * Y_field%pr(s)%cells(ii,jj,kk)
				rho_over_w_array(n) = rho_over_w_array(n) + &
					scalar_array(n) / this%thermo%thermo_ptr%molar_masses(s)
			end do
			if (s /= gamma) then
				sum_not_gamma = sum_not_gamma + &
					charm_face_value(scalar_array,velocity) / this%thermo%thermo_ptr%molar_masses(s)
			end if
		end do

		rho_over_w_face = charm_face_value(rho_over_w_array,velocity)
		phi = this%thermo%thermo_ptr%molar_masses(gamma) * &
			(rho_over_w_face - sum_not_gamma)
	end function

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
                        if (dimensions == 1) then
                            if (neighbours_distance(n) == dimensions) then
                                neighbours_coeffs(n) = 1.0_dp
                            else
                                neighbours_coeffs(n) = 0.0_dp
                            end if
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
                            if (neighbours_distance(n) == dimensions) then
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                                    !# Neumann boundary
                                    neighbours_coeffs(n) = (2.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)))
                                else
                                    !# Dirichlet boundary
                                    neighbours_coeffs(n) = (2.0_dp + neighbours_bound(nn(1)) + neighbours_bound(nn(2)))
                                end if
                            else if (neighbours_distance(n) == dimensions - 1) then
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                                    if (neighbours_distance(nn(1)) == dimensions) then
                                        !# Neumann boundary
                                        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(1)))
                                    else
                                        neighbours_coeffs(n) = (1.0_dp - neighbours_bound(nn(2)))
                                    end if
                                else
                                    !# Dirichlet boundary
                                    if (neighbours_distance(nn(1)) == dimensions) then
                                        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(1)))
                                    else
                                        neighbours_coeffs(n) = (1.0_dp + neighbours_bound(nn(2)))
                                    end if
                                end if
                            else
                                if ((neighbours_bound(nn(1)) == 1).or.(neighbours_bound(nn(2)) == 1)) then
                                    !# Neumann boundary
                                    neighbours_coeffs(n) = 0.0_dp
                                else
                                    !# Dirichlet boundary
                                    neighbours_coeffs(n) = 0.0_dp
                                end if
                            end if
                            
                            
                            neighbours_coeffs(n) = 1.0_dp / 4.0_dp * neighbours_coeffs(n)
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
		
		if (time > 1e-03) then
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
		delay	= 0.0e-03_dp
		duration = 5000.0e-06_dp
		
		associate (rho	=> this%rho%s_ptr)
			if ((time <= delay + duration)) then 
				iteration = iteration + 1
				do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
				do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
				do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
                    if (( (j - 5)**2 + (i - 200)**2 <= 5)) then
						rho%cells(i,j,:)	= 0.1_dp * this%rho_0 !- this%rho_0*18.5_dp/20.0_dp! *(time-delay)/relax
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
			if((T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1) > max_grad_T) then
				max_grad_T = (T%cells(i+1,1,1) - T%cells(i-1,1,1))/cell_size(1)
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
 
        name_string = 'VARIABLES='
		name_string = trim(name_string) // 'x'
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
