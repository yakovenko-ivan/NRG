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
	
	use solver_options_class
	
	implicit none
	
#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: fds_solver, fds_solver_c

	type(field_scalar_cons)	,target	:: p_dyn	,p_stat		,p_stat_old	,dp_stat_dt	,p_int, dT_dt, T_int, rho_int, rho_old	
	type(field_scalar_cons)	,target	:: div_v	,div_v_int	,ddiv_v_dt	,H	, H_old
	type(field_scalar_cons)	,target	:: E_f_int
	type(field_vector_cons)	,target	:: v_int, Y_int, Y_old
	type(field_scalar_flow)	,target	:: F_a, F_b
	type(field_vector_flow)	,target	:: v_f, v_f_old	

	real(dkind)	,dimension(:)	,allocatable	:: concs
	!$omp threadprivate(concs)

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

		type(lagrangian_droplets_solver), dimension(:)	    ,allocatable	:: droplets_solver		
		
		type(computational_domain)					:: domain
		type(thermophysical_properties_pointer)		:: thermo		
		type(chemical_properties_pointer)			:: chem
		type(computational_mesh_pointer)			:: mesh
		type(boundary_conditions_pointer)			:: boundary

		type(field_scalar_cons_pointer)	:: rho		, rho_int		, rho_old		, T				, T_int			, dT_dt			, p				, p_int			, v_s			, mol_mix_conc
		type(field_scalar_cons_pointer)	:: E_f		, E_f_prod_chem	, E_f_prod_heat	, E_f_prod_gd	, E_f_prod_visc	, E_f_prod_diff	, E_f_int		, h_s			, gamma
		type(field_scalar_cons_pointer)	:: p_stat	, p_stat_old	, dp_stat_dt	, p_dyn			, div_v			, div_v_int		, ddiv_v_dt		, H				, H_old			
		type(field_scalar_cons_pointer)	:: nu		, kappa
		type(field_scalar_flow_pointer)	:: F_a		, F_b
		
		type(field_vector_cons_pointer)	:: v		, v_prod_gd		, v_prod_visc	, v_prod_source	, v_int	
		type(field_vector_cons_pointer)	:: Y		, Y_prod_diff	, Y_prod_chem	, Y_int			, Y_old
		type(field_vector_cons_pointer)	:: D
		type(field_vector_flow_pointer)	:: v_f		, v_f_old
		
		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	::  rho_prod_droplets, E_f_prod_droplets
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	::  Y_prod_droplets		
		type(field_vector_flow_pointer)	,dimension(:)	,allocatable	::  v_prod_droplets
		
		real(dkind)	,dimension(:,:,:,:)	,allocatable	:: vorticity, grad_F_a, grad_F_b
		real(dkind)	,dimension(:,:,:)	,allocatable	:: p_old
		
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
		procedure	,private	:: spectral_radii_selection
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
		integer				:: particles_phase_counter, droplets_phase_counter		
		
		character(len=40)	:: var_name
		
		integer				:: bound_number, plus, sign
		character(len=20)	:: boundary_type_name

		real(dkind)	,dimension(3)	:: cell_size
		real(dkind)					:: x, y
		integer	:: dimensions
		integer	:: i,j,k,dim,dim1
		
		constructor%diffusion_flag		= problem_solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= problem_solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= problem_solver_options%get_heat_transfer_flag()
		constructor%reactive_flag		= problem_solver_options%get_chemical_reaction_flag()
		constructor%hydrodynamics_flag	= problem_solver_options%get_hydrodynamics_flag()
		constructor%courant_fraction	= problem_solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= problem_solver_options%get_CFL_condition_flag()
		constructor%sources_flag		= .false.
		
		constructor%additional_droplets_phases_number	= problem_solver_options%get_additional_droplets_phases_number()		
		
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
	
		call manager%create_scalar_field(dT_dt,'temperature_change','dT_dt')
		constructor%dT_dt%s_ptr						=> dT_dt
		
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
		
		if(constructor%viscosity_flag) then
			constructor%visc_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity')
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
				constructor%droplets_solver(droplets_phase_counter)	= lagrangian_droplets_solver_c(manager, droplets_params, droplets_phase_counter)
				write(var_name,'(A,I2.2)') 'energy_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%E_f_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%rho_prod_droplets(droplets_phase_counter)%s_ptr	=> scal_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_droplets', droplets_phase_counter
				call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,var_name)
				constructor%v_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_f_ptr%v_ptr	
				write(var_name,'(A,I2.2)') 'concentration_production_droplets', droplets_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_ptr,vect_ptr,tens_ptr,var_name)
				constructor%Y_prod_droplets(droplets_phase_counter)%v_ptr		=> vect_ptr%v_ptr                
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
											
		cons_utter_loop	= manager%domain%get_local_utter_cells_bounds()	
		cons_inner_loop = manager%domain%get_local_inner_cells_bounds()	

		problem_data_io				= data_io_c(manager,calculation_time)									
		
		call problem_data_io%input_all_data()
		
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
		end if		
		dimensions		= manager%domain%get_domain_dimensions()

		constructor%all_Neumann_flag = .true.
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(constructor%boundary%bc_ptr%bc_markers(i,j,k) == 0) then
				do dim = 1,dimensions	
					do plus = 1,2
						sign			= (-1)**plus
						bound_number	= constructor%boundary%bc_ptr%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
						if( bound_number /= 0 ) then
							boundary_type_name = constructor%boundary%bc_ptr%boundary_types(bound_number)%get_type_name()
							if (boundary_type_name == 'outlet') then
								constructor%all_Neumann_flag = .false.
							end if
						end if
					end do
				end do
			end if
		end do
		end do
		end do		


		cell_size						= constructor%mesh%mesh_ptr%get_cell_edges_length()
		
		constructor%time				= calculation_time
		constructor%initial_time_step	= 2.5e-08_dkind ! problem_solver_options%get_initial_time_step()
		constructor%time_step			= constructor%initial_time_step

		constructor%rho_0				= constructor%rho%s_ptr%cells(1,1,1)
	
		species_number = manager%chemistry%chem_ptr%species_number
		
		!$omp parallel
		allocate(concs(species_number))
		concs = 0.0_dkind
		!$omp end parallel
		
	end function

	subroutine solve_problem(this,iteration)
		class(fds_solver)	,intent(inout)	:: this
		integer				,intent(in)		:: iteration
		
		integer	:: droplets_phase_counter
		
		integer	:: specie



		this%time = this%time + this%time_step		

		if (this%reactive_flag)		call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
		if (this%viscosity_flag)	call this%visc_solver%solve_viscosity(this%time_step)
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)				
		
		if(this%additional_droplets_phases_number /= 0) then
			do droplets_phase_counter = 1, this%additional_droplets_phases_number
				call this%droplets_solver(droplets_phase_counter)%droplets_solve(this%time_step)		
			end do		
		end if  		
		
		call this%calculate_interm_Y_predictor(this%time_step)
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.true.)		
		call this%apply_boundary_conditions(this%time_step,predictor=.true.)		
		call this%calculate_divergence_v		(this%time_step,predictor=.true.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.true.)
		call this%calculate_velocity			(this%time_step,predictor=.true.)
		
		call this%calculate_interm_Y_corrector(this%time_step)
		call this%state_eq%apply_state_equation_low_mach_fds(this%time_step,predictor=.false.)
		call this%apply_boundary_conditions(this%time_step,predictor=.false.)

		call this%calculate_divergence_v		(this%time_step,predictor=.false.)
		call this%calculate_pressure_poisson	(this%time_step,predictor=.false.)
		call this%calculate_velocity			(this%time_step,predictor=.false.)		

		call this%calculate_time_step()

		!call this%state_eq%check_conservation_laws()

	end subroutine

	subroutine calculate_interm_Y_predictor(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		
		real(dkind)	:: B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dkind)	,dimension(:)	,allocatable	:: Y_rho_int
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
		
		allocate(Y_rho_int(species_number))
		
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
					
			!$omp parallel default(none)  private(i,j,k,dim,spec,spec_summ,Y_rho_int,rhs,flux_right,flux_left,lame_coeffs) , &
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
							! x -> z, y -> r
							lame_coeffs(2,1)	=  mesh%mesh(2,i,j,k) - 0.5_dkind*cell_size(1)			
							lame_coeffs(2,2)	=  mesh%mesh(2,i,j,k)
							lame_coeffs(2,3)	=  mesh%mesh(2,i,j,k) + 0.5_dkind*cell_size(1)	
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
						
						if (this%diffusion_flag)	Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_diff%pr(spec)%cells(i,j,k) * time_step
						if (this%reactive_flag)		Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_chem%pr(spec)%cells(i,j,k) * time_step						
						
						if (this%additional_droplets_phases_number /= 0) then
							do droplets_phase_counter = 1, this%additional_droplets_phases_number
								Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_droplets(droplets_phase_counter)%v_ptr%pr(spec)%cells(i,j,k)
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
	
		real(dkind)	:: flux_left, flux_right, specie_enthalpy, mixture_cp, mixture_cp_dT, dp_dt, D_sum, P_sum, U_sum, div_sum
		
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
					dT_dt			=> this%dT_dt%s_ptr			, &
					rho				=> this%rho%s_ptr			, &
					p_stat			=> this%p_stat%s_ptr		, &
					gamma			=> this%gamma%s_ptr			, &
					h_s				=> this%h_s%s_ptr			, &
					mol_mix_conc	=> this%mol_mix_conc%s_ptr	, &
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
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

			D_sum = 0.0_dkind
			P_sum = 0.0_dkind		
			U_sum = 0.0_dkind
			
			div_sum = 0.0_dkind
			
			time = time + 0.5_dkind*time_step
			
			if(time <= 2.0e-03_dkind) then
				energy_source = 1.0e08_dkind
			else
				energy_source = 0.0_dkind
			end if
			
		!	print *, '2', div_v_int%cells(50,20,1),T%cells(50,20,1),D_sum,P_sum
			
			!$omp parallel default(none)  private(flux_right,flux_left,i,j,k,dim,spec,mixture_cp,mixture_cp_dT,specie_enthalpy,plus,sign,bound_number,cell_volume,cell_surface_area,lame_coeffs) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(D_sum,P_sum,U_sum,div_v_int,p_stat,dp_stat_dt,v_f,h_s,E_f_prod_visc,E_f_prod_heat,E_f_prod_diff,E_f_prod_chem,Y_prod_diff,Y_prod_chem,T,dT_dt,Y,rho,mol_mix_conc,energy_source,bc,cons_inner_loop,species_number,dimensions,cell_size,coordinate_system,mesh,iter)					
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
					
					if (this%viscosity_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_visc%cells(i,j,k)
					if (this%heat_trans_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_heat%cells(i,j,k)
					if (this%diffusion_flag)	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_diff%cells(i,j,k)
					if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) +  E_f_prod_chem%cells(i,j,k)
					

					concs = 0.0_dkind
					mixture_cp = 0.0
					do spec = 1,species_number
						if (this%thermo%thermo_ptr%molar_masses(spec) /= 0.0_dkind) then
							mixture_cp		= mixture_cp + this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i,j,k),spec)	*Y%pr(spec)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(spec)
							concs(spec)		= Y%pr(spec)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(spec) * mol_mix_conc%cells(i,j,k)
						end if
					end do				
					
					if(mesh%mesh(1,i,j,k) <= mesh%mesh(1,1,j,k) + 4.0e-04_dkind) then
				!		div_v_int%cells(i,j,k)	= div_v_int%cells(i,j,k) + energy_source !* 1.0e+10 * this%time_step *  rho%cells(i,j,k) !* cell_size(1) ! * 4.0_dkind * Pi * mesh%mesh(1,i,j,k) * mesh%mesh(1,i,j,k)
					end if
						
					if (this%additional_droplets_phases_number /= 0) then
						do droplets_phase_counter = 1, this%additional_droplets_phases_number
							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + E_f_prod_droplets(droplets_phase_counter)%s_ptr%cells(i,j,k)
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
				
					div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) / (rho%cells(i,j,k) * h_s%cells(i,j,k))
					mixture_cp		= this%thermo%thermo_ptr%calculate_mixture_cp(T%cells(i,j,k), concs)
					mixture_cp_dT	= this%thermo%thermo_ptr%calculate_mixture_cp_dT(T%cells(i,j,k),concs)
					
					div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) - 1.0_dkind/mixture_cp*mixture_cp_dT*dT_dt%cells(i,j,k)	
					
					do spec = 1, species_number
					
						specie_enthalpy = (this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i,j,k),spec))*T%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(spec)
								
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

							div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) *  h_s%cells(i,j,k))) * &
																			  (	-  (	v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* lame_coeffs(dim,3) * (flux_right	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)	&
																					 -	v_f%pr(dim)%cells(dim,i,j,k)									* lame_coeffs(dim,1) * (flux_left	-  rho%cells(i,j,k) * Y%pr(spec)%cells(i,j,k)) / cell_size(1)) / lame_coeffs(dim,2))
							continue													 
						end do							
						
						
						if (this%diffusion_flag) 	div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) * h_s%cells(i,j,k))) * &
																		  (	Y_prod_diff%pr(spec)%cells(i,j,k))

						if (this%reactive_flag)		div_v_int%cells(i,j,k) = div_v_int%cells(i,j,k) + (	mol_mix_conc%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec) / rho%cells(i,j,k) - specie_enthalpy / (rho%cells(i,j,k) *  h_s%cells(i,j,k))) * &
																		  (	Y_prod_chem%pr(spec)%cells(i,j,k))
					end do
					
					D_sum = D_sum + div_v_int%cells(i,j,k) * cell_volume
					P_sum = P_sum + (1.0_dkind / p_stat%cells(i,j,k) - 1.0_dkind / (rho%cells(i,j,k) * h_s%cells(i,j,k)))* cell_volume 
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
					
					mixture_cp = 0.0
					do spec = 1, species_number
						mixture_cp		= mixture_cp + this%thermo%thermo_ptr%calculate_specie_cp(T%cells(i,j,k),spec)*Y%pr(spec)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(spec)
					end do


					dp_stat_dt%cells(i,j,k) = (D_sum - U_sum)/ P_sum
						
					if (this%all_Neumann_flag) then	
						div_v_int%cells(i,j,k)  = div_v_int%cells(i,j,k) - (1.0_dkind / p_stat%cells(i,j,k) - 1.0_dkind / (rho%cells(i,j,k) * h_s%cells(i,j,k)))* dp_stat_dt%cells(i,j,k)
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
		
		real(dkind)	:: H_center, H_left, H_right, F_a_left, F_a_right, F_b_left, F_b_right
		real(dkind)	:: H_residual, H_max, H_max_old, H_average, residual, a_norm_init, a_norm, a_norm_prev, H_summ, sum_ddiv_v_dt, grad_F_a_summ, grad_F_b_summ
		real(dkind)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		integer		:: r_i, r_j, r_k

		real(dkind)	,dimension(3)	:: cell_size		
		real(dkind)	:: time_step_adj, beta, beta_old, spectral_radii_jacobi, edges_number, cells_number
		real(dkind), save	:: time = 0.0_dkind
		
		integer	:: dimensions, iterations, iterations_number
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: droplets_phase_counter
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus
		
		integer			:: poisson_iteration , pressure_iteration , overall_poisson_iteration
		
		logical			:: converged = .false.
		logical			:: pressure_converged = .false.
		
		dimensions		= this%domain%get_domain_dimensions()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		
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
					vorticity		=> this%vorticity			, &
					grad_F_a		=> this%grad_F_a			, &
					grad_F_b		=> this%grad_F_b			, &
					v_f				=> this%v_f%v_ptr			, &
					v_f_old			=> this%v_f_old%v_ptr		, &
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					v_prod_droplets	=> this%v_prod_droplets		, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

		
			ddiv_v_dt%cells	= 0.0_dkind	
			sum_ddiv_v_dt	= 0.0_dkind
			vorticity		= 0.0_dkind	
			F_a%cells		= 0.0_dkind		
			F_b%cells		= 0.0_dkind		
			cells_number	= 0.0_dkind
			
			grad_F_a_summ	= 0.0_dkind
			
			!beta = 2.0_dkind/(1.0_dkind + sin(pi*cell_size(1)))
			!beta = 2.0_dkind - 2.0_dkind * pi*cell_size(1)
			!beta = 1.99927_dkind
			beta = 1.0
			!beta = 2.0_dkind/(1.0_dkind + sin(pi*cell_size(1)**2*cons_utter_loop(2,1)))
			
			spectral_radii_jacobi = 0.0_dkind
			do dim = 1, dimensions
				spectral_radii_jacobi = spectral_radii_jacobi + cell_size(1)**2 * cos(Pi/cons_utter_loop(dim,2)) / (dimensions * cell_size(1)**2)
			end do
			
			spectral_radii_jacobi =	    0.9999913706818E+00_dkind

			!$omp parallel default(none)  private(i,j,k,dim,dim1,dim2,loop,lame_coeffs,farfield_velocity) , &
			!$omp& firstprivate(this) , &
			!$omp& shared(ddiv_v_dt,sum_ddiv_v_dt,div_v_int,vorticity,v_f,v_f_old,F_a,grad_F_a, grad_F_a_summ, rho_old,rho_int,v_prod_visc,bc,cons_inner_loop,cons_utter_loop,flow_inner_loop,dimensions,predictor,cell_size,cells_number,coordinate_system,mesh,time_step)					
			!$omp do collapse(3) schedule(guided)	reduction(+:sum_ddiv_v_dt) reduction(+:cells_number)
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
							ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k)  - 0.5_dkind * ( (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / cell_size(1)	&
																							+ (v_f_old%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) * lame_coeffs(dim,3) - v_f_old%pr(dim)%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2) / cell_size(1)) / (0.5_dkind * time_step)
						end do
					end if
					
					cells_number  = cells_number + 1.0_dkind
					sum_ddiv_v_dt = sum_ddiv_v_dt + ddiv_v_dt%cells(i,j,k)
				!	ddiv_v_dt%cells(i,j,k)	= - 4.0_dkind*exp(-8.0**2*(i*1e-05 - 0.25)**2)
					
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
					ddiv_v_dt%cells(i,j,k)	= ddiv_v_dt%cells(i,j,k) - sum_ddiv_v_dt/cells_number
				end if
			end do
			end do
			end do		
			!$omp end do			
			end if
			
			!$omp do collapse(3) schedule(guided)	
			do k = flow_inner_loop(3,1),flow_inner_loop(3,2)
			do j = flow_inner_loop(2,1),flow_inner_loop(2,2)
			do i = flow_inner_loop(1,1),flow_inner_loop(1,2)
				do dim =  1, 3
				do dim1 = 1, dimensions
				do dim2 = 1, dimensions
					if((bc%bc_markers(i,j,k) == 0).and.(bc%bc_markers(i-I_m(dim1,1),j-I_m(dim1,2),k-I_m(dim1,3)) == 0).and.(bc%bc_markers(i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3)) == 0)) then 
						if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
							if(((dim1-dim) == 1).or.((dim1-dim) == -2)) then
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) - (v_f%pr(dim1)%cells(dim1,i,j,k) - v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))) / cell_size(1)
							else
								vorticity(dim,i,j,k) = vorticity(dim,i,j,k) + (v_f%pr(dim1)%cells(dim1,i,j,k) - v_f%pr(dim1)%cells(dim1,i-I_m(dim2,1),j-I_m(dim2,2),k-I_m(dim2,3))) / cell_size(1)
							end if
						end if
					end if
				end do
				end do
				end do
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
							F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(this%rho_0 - rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* g(dim))
						
							if (this%viscosity_flag)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(0.5_dkind*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))

							if (this%additional_droplets_phases_number /= 0) then
								do droplets_phase_counter = 1, this%additional_droplets_phases_number
									F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)
								end do		
							end if
							
							do dim1 = 1, 3
							do dim2 = 1, dimensions
								if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
								if(((dim-dim1) == 1).or.((dim-dim1) == -2)) then							
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								else
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								end if
								end if
							end do
							end do

						else
							F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_int%cells(i,j,k))) *(this%rho_0 - rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))* g(dim))

							if (this%viscosity_flag)	F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - (1.0_dkind/(0.5_dkind*(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + rho_old%cells(i,j,k))) *(0.5_dkind*(v_prod_visc%pr(dim)%cells(i,j,k) + v_prod_visc%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))))
							
							if (this%additional_droplets_phases_number /= 0) then
								do droplets_phase_counter = 1, this%additional_droplets_phases_number
									F_a%cells(dim,i,j,k) = F_a%cells(dim,i,j,k) + v_prod_droplets(droplets_phase_counter)%v_ptr%pr(dim)%cells(dim,i,j,k)
								end do
							end if							
							
							do dim1 = 1, 3
							do dim2 = 1, dimensions
								if ((dim1 /= dim).and.(dim2 /= dim1).and.(dim2 /= dim)) then
								if(((dim-dim1) == 1).or.((dim-dim1) == -2)) then							
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) - 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))
								else
									F_a%cells(dim,i,j,k)	=  F_a%cells(dim,i,j,k) + 0.25_dkind * (	vorticity(dim1,i,j,k)										* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_f%pr(dim2)%cells(dim2,i,j,k)) + &
																										vorticity(dim1,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))	* (v_f%pr(dim2)%cells(dim2,i-I_m(dim,1)+I_m(dim2,1),j-I_m(dim,2)+I_m(dim2,2),k-I_m(dim,3)+I_m(dim2,3)) + v_f%pr(dim2)%cells(dim2,i+I_m(dim2,1),j+I_m(dim2,2),k+I_m(dim2,3))))

								end if
								end if

							end do
							end do
						
						end if
					end if
				end do
				end do
				end do
				!$omp end do
			end do

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
					grad_F_a(dim,i,j,k)  = grad_F_a(dim,i,j,k)  - grad_F_a_summ/cells_number
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
				
				!$omp parallel default(none)  private(i,j,k,dim,loop,plus,sign,bound_number,boundary_type_name,residual,lame_coeffs,farfield_velocity) , &
				!$omp& firstprivate(this) , &
				!$omp& shared(F_a,F_b,grad_F_b,grad_F_b_summ,v_f,v_f_old,p_dyn,H,ddiv_v_dt,rho_old,rho_int,bc,predictor,cons_utter_loop,cons_inner_loop,dimensions,cell_size,coordinate_system,mesh,a_norm_init,time_step,cells_number)
				
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
								F_b%cells(dim,i,j,k)=	- (	p_dyn%cells(i,j,k)	*rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_old%cells(i,j,k))		&
														/	(rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_old%cells(i,j,k))	&
														*	(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	/ cell_size(1)
							else
								F_b%cells(dim,i,j,k)=	- (	p_dyn%cells(i,j,k)	*rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		&
														+	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	*rho_int%cells(i,j,k))		&
														/	(rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ rho_int%cells(i,j,k))	&
														*	(1.0_dkind/rho_int%cells(i,j,k)	- 1.0_dkind/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	/ cell_size(1)
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
						grad_F_b(dim,i,j,k)  = grad_F_b(dim,i,j,k)  - grad_F_b_summ/cells_number
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
					
						residual	= 0.0_dkind
						
						residual	= residual + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

						do dim = 1, dimensions
							residual	= residual +	(H%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
							residual	= residual -	(H%cells(i,j,k) - H%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)
						
							residual	= residual +	cell_size(1)*grad_F_a(dim,i,j,k)

							residual	= residual +	cell_size(1)*grad_F_b(dim,i,j,k)

							do plus = 1,2
								sign			= (-1)**plus
								bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
								if( bound_number /= 0 ) then
									boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
									select case(boundary_type_name)
										case('wall')
											farfield_velocity = 0.0_dkind
											if(predictor) then
												residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																						!- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))* cell_size(1)/time_step	
																																							
											else
												residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																						!- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																																						!		+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))* cell_size(1)/(0.5_dkind*time_step)	
											end if													
										case('inlet')
											farfield_velocity = this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
											if(predictor) then
												residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &	
																																						- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
											else
												residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																								+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &
																																						- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																																								+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))* cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
											end if													
									end select
								end if
							end do
						end do

						a_norm_init = a_norm_init + abs(residual)
						
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
				beta				= 1.0_dkind
				!	call this%spectral_radii_selection(time_step,predictor,a_norm_init,spectral_radii_jacobi)				
				!	pause	
				!$omp parallel default(none)  private(i,j,k,dim,residual,plus,sign,bound_number,boundary_type_name,lame_coeffs,farfield_velocity) , &
				!$omp& firstprivate(this) , &
				!$omp& shared(H_summ,H,a_norm,a_norm_init,a_norm_prev,H_max,H_max_old,H_average,F_a,F_b,p_dyn,rho_old,rho_int,beta,ddiv_v_dt,v_f,v_f_old,cons_inner_loop,cons_utter_loop,bc,dimensions,cell_size,coordinate_system,mesh,converged,predictor,time_step,time,poisson_iteration,spectral_radii_jacobi)				
				do while ((.not.converged).and.(poisson_iteration <= 50000))
				
					!$omp barrier
				
					!$omp master
					H_max	= 0.0_dkind
					a_norm	= 0.0_dkind
					converged = .true.
					!$omp end master
					
					!$omp barrier

				!# Parallel successive overrelaxation
					!$omp do collapse(3) schedule(guided) reduction(+:a_norm)
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						if(bc%bc_markers(i,j,k) == 0) then
							if((mod(i+j+k,2) == 1)) then	
			  
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

								residual = 0.0_dkind
								
								residual = residual + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)								

								do dim = 1, dimensions
								
									residual	= residual +	(H%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
									residual	= residual -	(H%cells(i,j,k) - H%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

									residual	= residual +	cell_size(1)*grad_F_a(dim,i,j,k) 

									residual	= residual +	cell_size(1)*grad_F_b(dim,i,j,k) 

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
													farfield_velocity = 0.0_dkind
													if(predictor) then
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																							
													else
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																								
																																								
													end if													
												case('inlet')
													farfield_velocity = this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													if(predictor) then
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &
																																								- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
													else
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &
																																								- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																																										+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))* cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
													end if													
											end select		
										end if
									end do
								end do								
								
								H%cells(i,j,k)	= H%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*residual
  
								a_norm = a_norm + abs(residual)
     
							end if
						end if
					end do
					end do
					end do				
					!$omp end do

					!$omp barrier

					!$omp master
					if(poisson_iteration == 0) then
						beta		= 1.0_dkind/(1.0_dkind - 0.5_dkind * spectral_radii_jacobi*spectral_radii_jacobi)
					else
						beta		= 1.0_dkind/(1.0_dkind - 0.25_dkind * spectral_radii_jacobi*spectral_radii_jacobi * beta)
					end if
					!$omp end master
					
					!$omp barrier					
				
					!$omp do collapse(3) schedule(guided) reduction(+:a_norm)			
					do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
					do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
					do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
						if(bc%bc_markers(i,j,k) == 0) then
							if((mod(i+j+k,2) == 0)) then	
	    
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
							
								residual = 0.0_dkind
								
								residual = residual + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)

								do dim = 1, dimensions
									residual	= residual +	(H%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
									residual	= residual -	(H%cells(i,j,k) - H%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)
									
									residual	= residual +	cell_size(1)*grad_F_a(dim,i,j,k)

									residual	= residual +	cell_size(1)*grad_F_b(dim,i,j,k) 

									do plus = 1,2
										sign			= (-1)**plus
										bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
										if( bound_number /= 0 ) then
											boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
											select case(boundary_type_name)
												case('wall')
													farfield_velocity = 0.0_dkind
													if(predictor) then
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)
																																							
																																							
													else
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
																																							
													end if													
												case('inlet')
													farfield_velocity = this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													if(predictor) then
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &	
																																								- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))* cell_size(1)/time_step) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
													else
														residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																										+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)   &	
																																								- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																																										+	v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))* cell_size(1)/(0.5_dkind*time_step)) * lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	
													end if													
											end select
										end if
									end do
								end do	
								
								H%cells(i,j,k)	= H%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*residual
								
								a_norm = a_norm + abs(residual)
								
							end if
						end if
					end do
					end do
					end do				
					!$omp end do

					!$omp barrier

					!$omp master

					if(poisson_iteration == 0) then
						beta		= 1.0_dkind/(1.0_dkind - 0.25_dkind * spectral_radii_jacobi*spectral_radii_jacobi * beta)
					end if

					if((a_norm_init < 1e-10).and.(poisson_iteration == 0)) then
					!	print *, a_norm_init, a_norm
						a_norm_init =  a_norm
					!	pause
					end if						
						
					if((mod(poisson_iteration,100) == 0).and.(a_norm_init > 1e-10)) then
						print *, poisson_iteration,  a_norm/a_norm_init, a_norm_init, a_norm, beta
					end if					
						
					if ((a_norm > 1.0e-02*a_norm_init).and.(converged)) converged = .false.	

					if (poisson_iteration < 100) converged = .false.
						
					if ((poisson_iteration > 500).and.(a_norm_prev < a_norm)) then
						converged = .true.
						print *, 'Not converged***************************************'
			!			pause
					end if
						
					if (a_norm < 1e-05) converged = .true.
						
					a_norm_prev = a_norm

					poisson_iteration	= poisson_iteration + 1

					!$omp end master
					
					!$omp barrier

				end do	
				!$omp end parallel
				
				print *, 'Poisson iteration:', poisson_iteration

				!print *,  H_max
				
				!pause

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
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  H%cells(i,j,k) - sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																		+	v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) / time_step) * cell_size(1)		
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  H%cells(i,j,k) - sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																		+	(v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) / time_step) * cell_size(1)
											end if
										case('outlet')
											if (sign == 1) then		!# Правая граница
												if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dkind) then		!# Выток, берутся значения слева
												
													if (predictor) then	
														H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/rho_old%cells(i,j,k)
													else
														H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/rho_int%cells(i,j,k)
													end if
													
													do dim1 = 1,dimensions
														H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) + 0.5_dkind*(	(0.5_dkind *( v_f%pr(dim1)%cells(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) + v_f%pr(dim1)%cells(dim1,i,j,k))) **2) 
													end do
												else																					!# Вток, берутся значения справа
													farfield_density		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_density()
													farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													farfield_pressure		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
							
													H%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))		=	p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/farfield_density + 0.5_dkind*(farfield_velocity **2)
											
												end if
											else					!# Левая граница
												if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dkind) then										!# Выток, берутся значения справа
										
													if (predictor) then	
														H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/rho_old%cells(i,j,k)
													else
														H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/rho_int%cells(i,j,k)
													end if
												
													do dim1 = 1,dimensions
														H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + 0.5_dkind*(	(0.5_dkind *( v_f%pr(dim1)%cells(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) + v_f%pr(dim1)%cells(dim1,i,j,k))) **2) 
													end do
												
												else																					!# Вток, берутся значения справа
													farfield_density		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_density()
													farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
													farfield_pressure		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
												
													H%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))		=	p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/farfield_density + 0.5_dkind*(farfield_velocity **2) 

												end if
											end if
												
										case('inlet')
     
											farfield_velocity = this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
												
											if (predictor) then
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = H%cells(i,j,k) - sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)	&	
																																- sign*(	farfield_velocity - v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)))* cell_size(1)/time_step
											else
												H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = H%cells(i,j,k) - sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))					&
																																		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) * cell_size(1)	&	
- sign*(	farfield_velocity - 0.5_dkind*(v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) &
																																		+ v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))))* cell_size(1)/(0.5_dkind*time_step)
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
				
				!$omp parallel default(none)  private(i,j,k,dim,H_residual,r_i,r_j,r_k,plus,sign,bound_number,boundary_type_name,farfield_density,farfield_velocity,farfield_pressure) , &
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
										*(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>10.0_dkind/cell_size(1)/cell_size(1)) then
									pressure_converged = .false.
									print *, 'Pressure error', abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_old%cells(i,j,k)	- 1.0_dkind/rho_old%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))
								end if
								end if
							else
								if (pressure_converged) then
								if (abs(((p_dyn%cells(i,j,k) - this%p_old(i,j,k)) - (p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - this%p_old(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))/cell_size(1) &
										*(1.0_dkind/rho_int%cells(i,j,k)	- 1.0_dkind/rho_int%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))/ cell_size(1))>10.0_dkind/cell_size(1)/cell_size(1)) then
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
	
			!$omp parallel default(none)  private(i,j,k,dim,vel_abs,plus,sign,bound_number,boundary_type_name,farfield_density,farfield_velocity,farfield_pressure) , &
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
										if (sign == 1) then		!# Правая граница
											if (v_f%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dkind) then		!# Выток, берутся значения слева
												p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.0_dkind 
											else																					!# Вток, берутся значения справа
												p_dyn%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.0_dkind
											end if
										else					!# Левая граница
											if (v_f%pr(dim)%cells(dim,i,j,k) <= 0.0_dkind) then										!# Выток, берутся значения справа
												p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	0.0_dkind 
											else																					!# Вток, берутся значения справа
												p_dyn%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	=	0.0_dkind
											end if
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
		
		character(len=20)		:: boundary_type_name
		
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus
		
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
										v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) = 0.0_dkind 
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
	
	subroutine calculate_interm_Y_corrector(this,time_step)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		
		real(dkind)	:: rhs, B_r, flux_left, flux_right, phi_loc, phi_up, r, s
		real(dkind)	,dimension(:)	,allocatable	:: Y_rho_int
		real(dkind)	:: spec_summ
		
		real(dkind), dimension (3,3):: lame_coeffs
		
		real(dkind)	,dimension(3)	:: cell_size		
		
		integer	:: dimensions, species_number
		integer	,dimension(3,2)	:: cons_inner_loop
		character(len=20)	:: coordinate_system
		
		integer	:: bound_number
		integer :: i,j,k,dim,spec,droplets_phase_counter
		
		dimensions		= this%domain%get_domain_dimensions()
		species_number	= this%chem%chem_ptr%species_number
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
				
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		allocate(Y_rho_int(species_number))
		
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

					
			!$omp parallel default(none)  private(i,j,k,dim,spec,spec_summ,Y_rho_int,rhs,flux_right,flux_left,lame_coeffs) , &
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
								Y_rho_int(spec) = Y_rho_int(spec) + Y_prod_droplets(droplets_phase_counter)%v_ptr%pr(spec)%cells(i,j,k)
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
										T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= farfield_temperature
										farfield_pressure												= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_pressure()
										
										if (predictor) then
											rho_int%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											rho_old%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											do spec = 1, species_number
												Y_int%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
												Y_old%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= Y%pr(spec)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
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
										
										farfield_velocity		= this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_velocity()
										
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
			s	= 1/r
			B_r	= s * (3*s + 1) / (s+1) / (s+1)	
		else
			B_r = 0.0_dkind
		end if

		!# SUPERBEE
		!if ( r /= 0.0_dkind) then
		!	B_r	= max(max(0.0_dkind,min(2.0_dkind*r,1.0_dkind)),min(r,2.0_dkind))
		!end if
		
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
				!	velocity_value	= velocity_value + abs(v%pr(dim)%cells(i,j,k))/(cell_size(1))		!#L1 norm
					velocity_value	= velocity_value + (v%pr(dim)%cells(i,j,k)/(cell_size(1)))**2		!#L2 norm
				end do

				if((velocity_value > 0.0_dkind).or.(abs(div_v_int%cells(i,j,k)) > 0.0_dkind)) then
				!	delta_t_interm = 1.0_dkind/(velocity_value + abs(div_v_int%cells(i,j,k)))		!# L1 norm
					delta_t_interm = 1.0_dkind/(sqrt(velocity_value) + abs(div_v_int%cells(i,j,k)))	!# L2 norm

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
		
		this%time_step	= min(2.5e-05_dkind, this%time_step)

	!	print *, this%time_step
		
		time = time + this%time_step
		
	!	this%time_step	= time_step
		
		end associate			

	end subroutine

	subroutine spectral_radii_selection(this,time_step,predictor,a_norm_init,spectral_radii_jacobi)
		class(fds_solver)	,intent(inout)	:: this
		real(dkind)			,intent(in)		:: time_step
		logical				,intent(in)		:: predictor
		real(dkind)			,intent(in)		:: a_norm_init
		real(dkind)			,intent(out)	:: spectral_radii_jacobi

		real(dkind)	:: H_residual, residual, a_norm, a_norm_prev
		real(dkind)	:: farfield_density, farfield_pressure, farfield_velocity
		
		real(dkind), dimension (3,3)	:: lame_coeffs
		character(len=20)				:: coordinate_system
		
		real(dkind)	,dimension(3)	:: cell_size		
		real(dkind)	:: time_step_adj, beta, beta_old, delta, x1, x2, xm, a, b
		real(dkind), save	:: time = 0.0_dkind
		real(dkind)	,dimension(3)	:: x
		integer		,dimension(3)	:: poisson_iteration
		
		
		integer	:: dimensions, iterations, iterations_number
		integer	,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
		integer	,dimension(3,2)	:: loop
		
		character(len=20)		:: boundary_type_name
		
		integer	:: iter, func_iter
		integer	:: sign, bound_number
		integer :: i,j,k,dim,dim1,dim2,spec,plus
		integer	:: io_unit
		
		integer			:: pressure_iteration , overall_poisson_iteration, poisson_iteration_1, poisson_iteration_2, poisson_iteration_m
		
		logical			:: converged = .false.
		logical			:: pressure_converged = .false.
		
		dimensions		= this%domain%get_domain_dimensions()
		
		coordinate_system	= this%domain%get_coordinate_system_name()
		
		cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		cons_utter_loop	= this%domain%get_local_utter_cells_bounds()
		
		cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()

		open(newunit = io_unit, file = 'spectral_radii.dat', status = 'unknown', form = 'formatted', position='append')
		
		associate (	p				=> this%p%s_ptr				, &
					p_stat			=> this%p_stat%s_ptr		, &
					p_dyn			=> this%p_dyn%s_ptr			, &
					rho_old			=> this%rho_old%s_ptr		, &
					rho_int			=> this%rho_int%s_ptr		, &
					ddiv_v_dt		=> this%ddiv_v_dt%s_ptr		, &
					F_a				=> this%F_a%s_ptr			, &
					F_b				=> this%F_b%s_ptr			, &
					grad_F_a		=> this%grad_F_a			, &
					grad_F_b		=> this%grad_F_b			, &					
					H				=> this%H%s_ptr				, &
					H_old			=> this%H_old%s_ptr			, &
					mesh			=> this%mesh%mesh_ptr		, &
					bc				=> this%boundary%bc_ptr)

			do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
			do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
			do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
				H_old%cells(i,j,k)		= H%cells(i,j,k)
			end do
			end do
			end do			

			print *, 'Calculating optimal spectral radii'
			
			a_norm_prev			= a_norm_init
			
			a = 0.9_dkind
			b = 0.999999_dkind

			!spectral_radii_jacobi = 0.0_dkind
			!do dim = 1, dimensions
			!	spectral_radii_jacobi = spectral_radii_jacobi + cell_size(1)**2 * cos(0.05_dkind*13.5_dkind*Pi/cons_utter_loop(dim,2)) / (dimensions * cell_size(1)**2)
			!end do			
			
			do iter = 1, 100

				x(1) = b - 0.25_dkind * (b - a)
				x(2) = a + 0.5_dkind  * (b - a)	
				x(3) = a + 0.25_dkind * (b - a)

				poisson_iteration = 0
				
				do func_iter = 1, 3
				
					spectral_radii_jacobi = x(func_iter)
				
					do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
					do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
					do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
						H%cells(i,j,k)		= H_old%cells(i,j,k)
					end do
					end do
					end do		

					converged = .false.

				!	beta		= 1.0_dkind/(1.0_dkind - 0.5_dkind * spectral_radii_jacobi*spectral_radii_jacobi)
					beta	= 1.0_dkind	
					
					do while ((.not.converged).and.(poisson_iteration(func_iter) <= 50000))
			
						a_norm	= 0.0_dkind
						converged = .true.
						

						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
								if((mod(i+j+k,2) == 1)) then	
			  
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

									residual = 0.0_dkind
								
									residual = residual + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)								
								
									do dim = 1, dimensions
								
										residual	= residual +	(H%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
										residual	= residual -	(H%cells(i,j,k) - H%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

										residual	= residual +	cell_size(1)*grad_F_a(dim,i,j,k) !(F_a%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_a%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)

										residual	= residual +	cell_size(1)*grad_F_b(dim,i,j,k) !(F_b%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_b%cells(dim,i,j,k) * lame_coeffs(dim,1)) /  lame_coeffs(dim,2)
								
										do plus = 1,2
											sign			= (-1)**plus
											bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											if( bound_number /= 0 ) then
												boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
												select case(boundary_type_name)
													case('wall')
														if(predictor) then
															residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	!- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																																!		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																																!		+	v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) / time_step) * cell_size(1)		
														else
															residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	!- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																																!		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																																!		+	(v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) / time_step) * cell_size(1)
														end if
												end select
											end if
										end do
									end do								
								
									H%cells(i,j,k)	= H%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*residual
     
									a_norm = a_norm + abs(residual)
     
								end if
							end if
						end do
						end do
						end do				

						if(poisson_iteration(func_iter) == 0) then
							beta		= 1.0_dkind/(1.0_dkind - 0.5_dkind * spectral_radii_jacobi*spectral_radii_jacobi)
						else
							beta		= 1.0_dkind/(1.0_dkind - 0.25_dkind * spectral_radii_jacobi*spectral_radii_jacobi * beta)
						end if
		
						do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
						do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
						do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
							if(bc%bc_markers(i,j,k) == 0) then
								if((mod(i+j+k,2) == 0)) then	
	    
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
							
									residual = 0.0_dkind
								
									residual = residual + cell_size(1)*cell_size(1)*ddiv_v_dt%cells(i,j,k)								

									do dim = 1, dimensions
										residual	= residual +	(H%cells(i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) - H%cells(i,j,k))* lame_coeffs(dim,3) / lame_coeffs(dim,2)
								
										residual	= residual -	(H%cells(i,j,k) - H%cells(i-i_m(dim,1),j-i_m(dim,2),k-i_m(dim,3)))* lame_coeffs(dim,1) / lame_coeffs(dim,2)

										residual	= residual +	cell_size(1)*grad_F_a(dim,i,j,k) !cell_size(1)*(F_a%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_a%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2)

										residual	= residual +	cell_size(1)*grad_F_b(dim,i,j,k) !cell_size(1)*(F_b%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * lame_coeffs(dim,3) - F_b%cells(dim,i,j,k) * lame_coeffs(dim,1)) / lame_coeffs(dim,2)

										do plus = 1,2
											sign			= (-1)**plus
											bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
											if( bound_number /= 0 ) then
												boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
												select case(boundary_type_name)
													case('wall')
														if(predictor) then
															residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	!- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																																		!		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																																		!		+	v_f%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3)) / time_step) * cell_size(1)		
														else
															residual = residual + (H%cells(i,j,k) - H%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)))* lame_coeffs(dim,2+sign) / lame_coeffs(dim,2)	!- sign*(	F_a%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))		&			
																																																		!		+	F_b%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))				&		
																																																		!		+	(v_f_old%pr(dim)%cells(dim,i+max(sign,0)*I_m(dim,1),j+max(sign,0)*I_m(dim,2),k+max(sign,0)*I_m(dim,3))) / time_step) * cell_size(1)
														end if
												end select
											end if
										end do
									end do	
								
									H%cells(i,j,k)	= H%cells(i,j,k) + 1.0_dkind/(2.0_dkind*dimensions)*beta*residual
								
									a_norm = a_norm + abs(residual)
								
								end if
							end if
						end do
						end do
						end do								

						if(poisson_iteration(func_iter) == 0) then
							beta		= 1.0_dkind/(1.0_dkind - 0.25_dkind * spectral_radii_jacobi*spectral_radii_jacobi * beta)
						end if
						
						if((mod(poisson_iteration(func_iter),100) == 0).and.(a_norm_init > 1e-10)) then
							print *, poisson_iteration(func_iter),  a_norm/a_norm_init, a_norm_init, a_norm, beta, spectral_radii_jacobi
						end if	
						
						if (poisson_iteration(func_iter) < 10) converged = .false.
						if ((poisson_iteration(func_iter) > 10000).and.(a_norm_prev < a_norm)) then
							print *, 'Not converged'
							converged = .true.
						end if
						
						if ((a_norm > 1.0e-01*a_norm_init).and.(converged)) converged = .false.	
						
						if (a_norm < 1e-02) converged = .true.

						poisson_iteration(func_iter) = poisson_iteration(func_iter) + 1

						if(func_iter>1) then
						if((poisson_iteration(func_iter) > maxval(poisson_iteration(:func_iter-1)))) then
							exit
						end if
						end if
						
					end do	
				
				end do	
					
				
				write(*,'(A,I2,A,E14.7,A,I4,A,E14.7,A,I4,A,E14.7,A,I4)') 'iter: ', iter, '	x1: ', x(3), '  ', poisson_iteration(3), ' x2: ', x(1), '  ', poisson_iteration(1), ' xm: ', x(2), '  ', poisson_iteration(2)						
				write(io_unit,'(E20.13,I6)')	x(MINLOC(poisson_iteration)), minval(poisson_iteration)
				
				
				if (poisson_iteration(3) < poisson_iteration(2)) then
					b = x(2)
					x(2) = x(3)
				else
					if (poisson_iteration(1) < poisson_iteration(2)) then
						a = x(2)
						x(2) = x(1)
					else
						a = x(3)
						b = x(1)
					end if
				end if
			end do
			
			close(io_unit)
			
			spectral_radii_jacobi = x(1)
			
			print *, 'Optimal spectral radii', spectral_radii_jacobi
			
			do k = cons_utter_loop(3,1),cons_utter_loop(3,2)
			do j = cons_utter_loop(2,1),cons_utter_loop(2,2)
			do i = cons_utter_loop(1,1),cons_utter_loop(1,2)
				H%cells(i,j,k)		= H_old%cells(i,j,k)
			end do
			end do
			end do				

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
	
