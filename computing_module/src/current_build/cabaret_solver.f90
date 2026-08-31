!> @file cabaret_solver_refactored_documented.f90
!! @brief Production CABARET solver for compressible multicomponent reactive flow.
!!
!! The implementation combines a conservative finite-volume balance step with
!! characteristic reconstruction of face variables.  Cell-centred conservative
!! variables and face-centred flow variables are advanced on a compact stencil.
!!
!! The production material-interface treatment reconstructs mass fractions and
!! temperature with the contact family, obtains pressure and normal velocity from
!! a thermally-perfect Riemann problem, and then closes face density through the
!! mixture equation of state.  Physical source rates are coupled consistently to
!! the predictor, corrector and source-shifted characteristic bounds.
!!
!! This file is a refactoring of the validated v5 solver.  Disabled research
!! branches, one-off diagnostics and unused storage have been removed without
!! changing the active numerical algorithm.

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
	use riemann_solver_class
	use thermophysical_properties_class	
	use chemical_properties_class

	use viscosity_solver_class
	use fickean_diffusion_solver_class
	use fourier_heat_transfer_solver_class
	use thermal_radiation_solver_class
	use chemical_kinetics_solver_class

	use dispersed_phase_solver_class, only: dispersed_phase_solver, &
		dispersed_phase_solver_c
    
	use mpi_communications_class

	use benchmarking   
	use solver_options_class	
	
	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: cabaret_solver, cabaret_solver_c

	! -------------------------------------------------------------------------
	! Production numerical constants
	! -------------------------------------------------------------------------
	! A face is treated as a material interface when adjacent mixture molar
	! masses differ by more than this ratio.
	real(dp), parameter :: material_contact_molar_mass_ratio = 1.01_dp

	! Exact-contact shortcut tolerances.  When p, u and T are uniform across a
	! composition jump, the exact star values are their arithmetic means.
	real(dp), parameter :: contact_equilibrium_p_rel_tol = 1.0e-8_dp
	real(dp), parameter :: contact_equilibrium_u_rel_tol = 1.0e-8_dp
	real(dp), parameter :: contact_equilibrium_T_rel_tol = 1.0e-8_dp

	! Initial face states use a thermally-perfect Riemann solve only where the
	! adjacent states are discontinuous according to these scale-free tests.
	real(dp), parameter :: initial_riemann_molar_mass_ratio = 1.01_dp
	real(dp), parameter :: initial_riemann_pressure_ratio = 1.05_dp
	real(dp), parameter :: initial_riemann_velocity_ratio = 0.05_dp
	real(dp), parameter :: initial_contact_velocity_eps = 10.0_dp*tiny(1.0_dp)

	! Material-contact reconstruction and admissible-state projection.
	real(dp), parameter :: contact_reconstruction_molar_mass_ratio = 1.000001_dp
	real(dp), parameter :: direct_psi_contact_molar_mass_ratio = 1.05_dp
	real(dp), parameter :: direct_psi_face_mixture_rel_tol = 2.0e-2_dp
	real(dp), parameter :: direct_psi_temperature_rel_tol = 3.0e-2_dp
	real(dp), parameter :: direct_psi_temperature_abs_tol = 50.0_dp

	! Final density invariant-domain bounds.  The absolute tolerance prevents
	! accidental activation for nearly vacuum states, while the relative part
	! leaves physical shock compression unaltered.
	real(dp), parameter :: face_density_guard_molar_mass_ratio = 1.01_dp
	real(dp), parameter :: face_density_guard_rel_tol = 0.10_dp
	real(dp), parameter :: face_density_guard_abs_tol = 1.0_dp

	type(field_scalar_flow)	,target	::	rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new, T_f_new
	type(field_vector_flow)	,target	::	Y_f_new, v_f_new
    
#ifdef OMP
	integer(kind=omp_lock_kind)	,dimension(:,:,:)	,allocatable	:: lock
#endif	
	
    type(timer)     :: cabaret_timer
    type(timer)     :: cabaret_gas_dynamics_timer
    type(timer)     :: cabaret_eos_timer
    type(timer)     :: cabaret_chemistry_timer
    type(timer)     :: cabaret_diffusion_timer
    type(timer)     :: cabaret_heattransfer_timer
    type(timer)     :: cabaret_radiation_timer
    type(timer)     :: cabaret_viscosity_timer

	type cabaret_solver
		logical :: diffusion_flag, viscosity_flag, heat_trans_flag, radiation_flag
		logical :: reactive_flag, CFL_condition_flag
		real(dp) :: courant_fraction
		real(dp) :: time, time_step
		real(dp)                    :: rho_0
		real(dp)    , dimension(3)  :: g
		integer			:: additional_particles_phases_number
        
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(diffusion_solver)				:: diff_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(thermal_radiation_solver)	:: radiation_solver
		type(viscosity_solver)				:: viscosity_solver	
		type(table_approximated_real_gas)	:: state_eq

		type(dispersed_phase_solver), dimension(:), allocatable :: particles_solver
		
        type(computational_domain)				:: domain
		type(mpi_communications)				:: mpi_support
		type(chemical_properties_pointer)		:: chem
		type(thermophysical_properties_pointer)	:: thermo
		type(computational_mesh_pointer)		:: mesh
		type(boundary_conditions_pointer)		:: boundary

		type(field_scalar_cons_pointer)	:: rho, T, p, v_s, gamma, E_f, mix_mol_mass
		type(field_scalar_flow_pointer)	:: rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new
		
		type(field_scalar_cons_pointer)	:: E_f_prod_chem, E_f_prod_heat, E_f_prod_rad, E_f_prod_diff, E_f_prod_visc

		type(field_vector_cons_pointer)	:: v, Y	, v_prod_visc
		type(field_vector_flow_pointer)	:: v_f_new, Y_f_new
		
		type(field_vector_cons_pointer)	:: Y_prod_chem, Y_prod_diff

		type(field_scalar_cons_pointer)	,dimension(:)	,allocatable	:: rho_prod_particles, E_f_prod_particles
		type(field_vector_cons_pointer)	,dimension(:)	,allocatable	:: Y_prod_particles
		type(field_vector_cons_pointer), dimension(:), allocatable :: v_prod_particles
        
		! Beginning-of-step conservative and thermodynamic state.
		real(dp), dimension(:,:,:), allocatable :: rho_old, p_old, E_f_old, v_s_old
		real(dp), dimension(:,:,:,:), allocatable :: v_old, Y_old

		! Effective conservative source rates evaluated once per full time step.
		real(dp), dimension(:,:,:), allocatable :: rho_src_old, rhoE_src_old
		real(dp), dimension(:,:,:,:), allocatable :: mom_src_old, rhoY_src_old

		! Old-time flux-divergence residuals used by the trapezoidal corrector.
		real(dp), dimension(:,:,:), allocatable :: rho_rhs_old, rhoE_rhs_old
		real(dp), dimension(:,:,:,:), allocatable :: mom_rhs_old, rhoY_rhs_old

		! Flow variables
		real(dp) ,dimension(:,:,:,:,:)	,allocatable    :: v_f, Y_f
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: rho_f, p_f, E_f_f, v_s_f
		! Quasi invariants
		real(dp) ,dimension(:,:,:,:,:)		,allocatable    :: r_inv_corr, q_inv_corr
        real(dp) ,dimension(:,:,:,:,:,:)		,allocatable    :: v_inv_corr

		! Contact-family thermodynamic coordinates at faces.
		! material_molar_volume_f=M(Y)/rho is the production density-closure
		! variable for strong multicomponent contacts.
		! material_sensible_energy_f is a dimensionless molar sensible
		! internal-energy coordinate used by the admissible pressure-equilibrium-preserving
		! material-contact correction.
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: material_temperature_f, material_molar_volume_f
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: material_sensible_energy_f
        
        
	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure	,private	:: apply_boundary_conditions_flow
		procedure				:: solve_problem
		procedure				:: calculate_time_step
		procedure				:: get_time_step
		procedure				:: get_time
		procedure				:: set_CFL_coefficient
		procedure	,private	:: cache_conservative_state
		procedure	,private	:: exchange_conservative_state
		procedure	,private	:: exchange_face_pressure_density
		procedure	,private	:: exchange_face_primitive_state
		procedure	,private	:: exchange_face_thermodynamic_state
		procedure	,private	:: zero_new_flow_state
		procedure	,private	:: normalize_face_mass_fractions
		procedure	,private	:: prepare_gas_dynamics_step
		procedure	,private	:: evaluate_effective_physical_source_rates
		procedure	,private	:: physical_source_acoustic_shifts
		procedure	,private	:: physical_source_contact_shifts
		procedure	,private	:: predict_conservative_half_step
		procedure	,private	:: update_cell_thermodynamics
		procedure	,private	:: reconstruct_acoustic_face_state
		procedure	,private	:: reconstruct_contact_family_face_state
		procedure	,private	:: finish_face_reconstruction
		procedure	,private	:: initialize_material_contact_faces
		procedure	,private	:: close_material_contact_density
		procedure	,private	:: apply_material_contact_riemann_closure
		procedure	,private	:: update_flow_thermodynamics
		procedure	,private	:: correct_conservative_full_step
		procedure	,private	:: finalize_gas_dynamics_step
		procedure	,private	:: advance_particle_phases
		procedure	,private	:: cache_flow_state_for_next_step
	end type

	interface	cabaret_solver_c
		module procedure	constructor
	end interface

contains
	!> Linearized contact-family quasi-invariant used by CABARET.
	pure real(dp) function contact_quasi_invariant(p_state, rho_state, c_state) result(inv)
		real(dp), intent(in) :: p_state, rho_state, c_state
		inv = p_state - c_state*c_state*rho_state
	end function contact_quasi_invariant

	!> Recover density from pressure and the linearized contact invariant.
	pure real(dp) function density_from_contact_quasi_invariant(p_state, inv, rho_ref, c_ref) result(rho_state)
		real(dp), intent(in) :: p_state, inv, rho_ref, c_ref
		rho_state = (p_state - inv)/(c_ref*c_ref)
	end function density_from_contact_quasi_invariant
	pure real(dp) function finite_volume_geometry_coefficient(nu, radius, dx) result(coeff)
		integer, intent(in) :: nu
		real(dp), intent(in) :: radius, dx
		real(dp) :: r_minus, r_plus

		coeff = 0.0_dp
		if (nu <= 1) return

		r_minus = radius - 0.5_dp*dx
		r_plus  = radius + 0.5_dp*dx

		! Same finite-volume geometric factor as in the conservative
		! predictor/corrector update. Invalid radial metrics are intentionally not
		! hidden by a numerical floor.
		coeff = 2.0_dp*real(nu - 1, dp)/(r_plus**(nu - 1) + r_minus**(nu - 1)) * &
			(r_plus**(nu - 1) - r_minus**(nu - 1))/(r_plus - r_minus)
	end function finite_volume_geometry_coefficient


	pure real(dp) function limit_quasi_invariant(value, left_value, half_value, right_value, source_shift, alpha, &
		& lower_clip, upper_clip) result(limited)
		real(dp), intent(in) :: value, left_value, half_value, right_value
		real(dp), intent(in) :: source_shift, alpha, lower_clip, upper_clip
		real(dp) :: min_inv, max_inv, width

		! Source-shifted maximum-principle interval following the
		! balance-characteristic form: the characteristic source term is not
		! interpreted as a spurious new extremum.
		min_inv = min(left_value, half_value, right_value) + source_shift
		max_inv = max(left_value, half_value, right_value) + source_shift
		width = abs(max_inv - min_inv)

		max_inv = max_inv + (-alpha)*width
		min_inv = min_inv - (-alpha)*width

		min_inv = max(min_inv, lower_clip)
		max_inv = min(max_inv, upper_clip)
		if (max_inv < min_inv) then
			max_inv = min_inv
		end if

		limited = min(max(value, min_inv), max_inv)
	end function limit_quasi_invariant


	!> Construct the solver, bind NRG fields, allocate work arrays and initialize face states.


	type(cabaret_solver)	function constructor(manager,problem_data_io)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io

		real(dp)								:: calculation_time
		
		type(field_scalar_cons_pointer)	:: scal_c_ptr
		type(field_vector_cons_pointer)	:: vect_c_ptr
		type(field_tensor_cons_pointer)	:: tens_c_ptr		

		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr
		type(field_tensor_flow_pointer)	:: tens_f_ptr	
        
        type(particles_phase)           :: particles_params
		integer							:: particles_phase_counter		
        
        character(len=40)       :: var_name
        
		real(dp)				:: spec_summ
		type(riemann_solver) :: initial_riemann
		real(dp)				:: molar_denom_left, molar_denom_right
		real(dp)				:: molar_mass_left, molar_mass_right, molar_mass_ratio
		real(dp)				:: rho_l, rho_r, p_l, p_r, v_l, v_r
		real(dp)				:: c_l, c_r, p_ratio, velocity_scale
		real(dp)				:: p_floor, rho_floor, u_face_init, u_contact_init
		real(dp)	,dimension(:)	,allocatable :: Y_left_riemann, Y_right_riemann, Y_face_riemann
		logical				:: use_left_contact_state, use_riemann_initial_face
		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds
		integer	,dimension(3,2)	:: cons_inner_loop, flow_inner_loop, loop
		
		integer					:: dimensions, species_number
		integer					:: i, j, k, dim, dim1, spec
		real(dp)	,dimension(3)	:: cell_size

		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		cons_inner_loop			= manager%domain%get_local_inner_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()

		species_number			= manager%chemistry%chem_ptr%species_number
		allocate(Y_left_riemann(species_number), Y_right_riemann(species_number), Y_face_riemann(species_number))
		
		cell_size				= manager%computational_mesh_pointer%mesh_ptr%get_cell_edges_length()
		
		constructor%diffusion_flag		= manager%solver_options%get_molecular_diffusion_flag()
		constructor%viscosity_flag		= manager%solver_options%get_viscosity_flag()
		constructor%heat_trans_flag		= manager%solver_options%get_heat_transfer_flag()
		constructor%radiation_flag		= manager%solver_options%get_thermal_radiation_flag()
		constructor%reactive_flag		= manager%solver_options%get_chemical_reaction_flag()
		constructor%courant_fraction	= manager%solver_options%get_CFL_condition_coefficient()
		constructor%CFL_condition_flag	= manager%solver_options%get_CFL_condition_flag()
        
        constructor%g                       = manager%solver_options%get_grav_acc()

        constructor%additional_particles_phases_number &
        	& = manager%solver_options%get_additional_particles_phases_number()
        
		constructor%domain				= manager%domain
		constructor%mpi_support			= manager%mpi_communications
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%boundary%bc_ptr		=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr
        
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'density')
		constructor%rho%s_ptr				=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'temperature')
		constructor%T%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure')
		constructor%p%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'full_energy')
		constructor%E_f%s_ptr				=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'internal_energy')
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_molar_mass')
		constructor%mix_mol_mass%s_ptr		=> scal_c_ptr%s_ptr		
		
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
		
		call manager%create_vector_field(Y_f_new,'specie_mass_fraction_flow'	,'Y_f_new',	'chemical')
		constructor%Y_f_new%v_ptr => Y_f_new		
		call manager%create_vector_field(v_f_new,'velocity_flow'					,'v_f_new',	'spatial')
		constructor%v_f_new%v_ptr => v_f_new	
		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity')
		constructor%v%v_ptr				=> vect_c_ptr%v_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_mass_fraction')
		constructor%Y%v_ptr				=> vect_c_ptr%v_ptr	
		
		if (constructor%reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_c_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_c_ptr%v_ptr
		end if

		if (constructor%diffusion_flag) then
			constructor%diff_solver			= diffusion_solver_c(manager, soret_enabled = manager%solver_options%get_soret_diffusion_flag())
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_diffusion')
			constructor%E_f_prod_diff%s_ptr			=> scal_c_ptr%s_ptr			
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_c_ptr%v_ptr
		end if

		if (constructor%heat_trans_flag .or. &
			constructor%additional_particles_phases_number > 0) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr			=> scal_c_ptr%s_ptr
		end if

		if (constructor%radiation_flag) then
			constructor%radiation_solver = thermal_radiation_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr, &
				'energy_production_radiation')
			constructor%E_f_prod_rad%s_ptr => scal_c_ptr%s_ptr
		end if
		if (constructor%viscosity_flag .or. &
			constructor%additional_particles_phases_number > 0) then
			constructor%viscosity_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_viscosity')
			constructor%E_f_prod_visc%s_ptr			=> scal_c_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_c_ptr%v_ptr
        end if			
		
        
		constructor%state_eq	=	table_approximated_real_gas_c(manager)        
        
		if(constructor%additional_particles_phases_number /= 0) then
			allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
			allocate(constructor%rho_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%Y_prod_particles(constructor%additional_particles_phases_number))
			do particles_phase_counter = 1, constructor%additional_particles_phases_number
				particles_params = manager%solver_options%get_particles_params(particles_phase_counter)
				constructor%particles_solver(particles_phase_counter) = dispersed_phase_solver_c( &
					manager, particles_params, particles_phase_counter)
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%rho_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter						
				call manager%get_cons_field_pointer_by_name( &
					scal_c_ptr, vect_c_ptr, tens_c_ptr, var_name)
				constructor%v_prod_particles(particles_phase_counter)%v_ptr => vect_c_ptr%v_ptr
				write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%Y_prod_particles(particles_phase_counter)%v_ptr		=> vect_c_ptr%v_ptr                
			end do		
		end if
	
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,tens_f_ptr,'adiabatic_index_flow')
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'adiabatic_index')
		constructor%gamma%s_ptr			=> scal_c_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_c_ptr%s_ptr		
		
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
								
									
									
									
		allocate(constructor%v_s_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
									
		
		allocate(constructor%v_old(		dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

									
		allocate(constructor%Y_old(		species_number					, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
								

		allocate(constructor%rho_src_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%rhoE_src_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%mom_src_old(	dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%rhoY_src_old(	species_number					, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%rho_rhs_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%rhoE_rhs_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%mom_rhs_old(	dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%rhoY_rhs_old(	species_number					, &
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

		allocate(constructor%material_temperature_f(	dimensions						, &
									flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
									flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
									flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		allocate(constructor%material_molar_volume_f(dimensions						, &
									flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
									flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
									flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		allocate(constructor%material_sensible_energy_f(dimensions			, &
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
        
		allocate(constructor%v_inv_corr(			2						, &	
											dimensions						, &
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
			constructor%material_temperature_f(:,:,:,:) = 0.0_dp
			constructor%material_molar_volume_f(:,:,:,:) = 0.0_dp
			constructor%material_sensible_energy_f(:,:,:,:) = 0.0_dp
			constructor%rho_src_old(:,:,:) = 0.0_dp
			constructor%rhoE_src_old(:,:,:) = 0.0_dp
			constructor%mom_src_old(:,:,:,:) = 0.0_dp
			constructor%rhoY_src_old(:,:,:,:) = 0.0_dp
			constructor%rho_rhs_old(:,:,:) = 0.0_dp
			constructor%rhoE_rhs_old(:,:,:) = 0.0_dp
			constructor%mom_rhs_old(:,:,:,:) = 0.0_dp
			constructor%rhoY_rhs_old(:,:,:,:) = 0.0_dp

			do dim = 1, dimensions		

				loop = flow_inner_loop

				do dim1 = 1,dimensions
					loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
				end do

				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2) 
				do i = loop(1,1),loop(1,2) 


							constructor%p_f(dim,i,j,k)		=	0.5_dp * (constructor%p%s_ptr%cells(i,j, &
								& k)	+ constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%rho_f(dim,i,j,k)	=	0.5_dp * (constructor%rho%s_ptr%cells(i,j, &
								& k)	+ constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%E_f_f(dim,i,j,k)	=	0.5_dp * (constructor%E_f%s_ptr%cells(i,j, &
								& k)	+ constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%v_s_f(dim,i,j,k)	=	0.5_dp * (constructor%v_s%s_ptr%cells(i,j, &
								& k)	+ constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							spec_summ = 0.0_dp
							do spec = 1, species_number
								constructor%Y_f(spec,dim,i,j,k) = 0.5_dp * (constructor%Y%v_ptr%pr(spec)%cells(i,j, &
									& k) + constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								spec_summ = spec_summ + max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp)
						end do

							do spec = 1,species_number
								constructor%Y_f(spec,dim,i,j,k) = max(constructor%Y_f(spec,dim,i,j,k), 0.0_dp) / spec_summ
							end do

							do dim1 = 1, dimensions
								constructor%v_f(dim1,dim,i,j,k) = 0.5_dp * (constructor%v%v_ptr%pr(dim1)%cells(i,j, &
									& k) + constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) )
							end do

							! If the initial face lies on a discontinuity, do not
							! initialize the CABARET flux variables by arithmetic averaging.
							! The first predictor step uses these old face states before any
							! characteristic reconstruction is available.  At a H2/air contact,
							! an averaged Y_f produces a nonphysical mixture face and can
							! generate a temperature spike immediately.  For such faces,
							! initialize the complete face state from a local 1-D Riemann
							! solution and select species/tangential velocity from the
							! material side indicated by the contact speed.
							if ((constructor%boundary%bc_ptr%bc_markers(i,j,k) == 0) .and. &
								(constructor%boundary%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

								molar_denom_left  = 0.0_dp
								molar_denom_right = 0.0_dp
								do spec = 1, species_number
								molar_denom_left = molar_denom_left + &
									constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) / &
									constructor%thermo%thermo_ptr%molar_masses(spec)
								molar_denom_right = molar_denom_right + &
									constructor%Y%v_ptr%pr(spec)%cells(i,j,k) / &
									constructor%thermo%thermo_ptr%molar_masses(spec)
								end do

								molar_mass_left  = 1.0_dp / molar_denom_left
								molar_mass_right = 1.0_dp / molar_denom_right
								molar_mass_ratio = max(molar_mass_left,molar_mass_right) / &
									min(molar_mass_left,molar_mass_right)

								rho_l   = constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								rho_r   = constructor%rho%s_ptr%cells(i,j,k)
								p_l     = constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								p_r     = constructor%p%s_ptr%cells(i,j,k)
								v_l     = constructor%v%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								v_r     = constructor%v%v_ptr%pr(dim)%cells(i,j,k)
								c_l     = constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								c_r     = constructor%v_s%s_ptr%cells(i,j,k)

								p_ratio = max(p_l,p_r) / min(p_l,p_r)
								velocity_scale = 0.5_dp*(abs(c_l) + abs(c_r))
								use_riemann_initial_face = (&
									 molar_mass_ratio > initial_riemann_molar_mass_ratio .or. &
									 p_ratio > initial_riemann_pressure_ratio .or. &
									 abs(v_l - v_r) > initial_riemann_velocity_ratio*velocity_scale)

								if (use_riemann_initial_face) then
									p_floor = 1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r))
									rho_floor = 1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r))
									do spec = 1, species_number
										Y_left_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
										Y_right_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
									end do
										call initial_riemann%set_thermally_perfect_parameters(constructor%thermo%thermo_ptr, &
											rho_l, rho_r, p_l, p_r, v_l, v_r, Y_left_riemann, Y_right_riemann, p_floor, rho_floor)
									call initial_riemann%solve()

									if (initial_riemann%get_success()) then
										u_face_init = initial_riemann%get_velocity()
										u_contact_init = initial_riemann%get_contact_velocity()
										use_left_contact_state = (u_contact_init >= 0.0_dp)
										call initial_riemann%get_mass_fractions(Y_face_riemann)

										constructor%p_f(dim,i,j,k)   = initial_riemann%get_pressure()
										constructor%rho_f(dim,i,j,k) = initial_riemann%get_density()
										constructor%v_f(dim,dim,i,j,k) = u_face_init
									else
										! Robust fallback: keep a single material side rather than an
										! averaged composition.  At an initially motionless diaphragm,
										! the pressure jump gives the first contact direction.
										u_face_init = 0.5_dp*(v_l + v_r)
										if (abs(u_face_init) > initial_contact_velocity_eps) then
											use_left_contact_state = (u_face_init >= 0.0_dp)
										else
											use_left_contact_state = (p_l >= p_r)
										end if
										if (use_left_contact_state) then
											constructor%p_f(dim,i,j,k) = p_l
											constructor%rho_f(dim,i,j,k) = rho_l
											constructor%v_f(dim,dim,i,j,k) = v_l
										else
											constructor%p_f(dim,i,j,k) = p_r
											constructor%rho_f(dim,i,j,k) = rho_r
											constructor%v_f(dim,dim,i,j,k) = v_r
										end if
									end if

									if (use_left_contact_state) then
										constructor%E_f_f(dim,i,j,k) = constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
										constructor%v_s_f(dim,i,j,k) = constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
										do spec = 1, species_number
											if (initial_riemann%get_success()) then
												constructor%Y_f(spec,dim,i,j,k) = Y_face_riemann(spec)
											else
												constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
													& 3))
											end if
										end do
										do dim1 = 1, dimensions
										if (dim1 /= dim) then
											constructor%v_f(dim1,dim,i,j,k) = constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
												& 3))
										end if
										end do
									else
										constructor%E_f_f(dim,i,j,k) = constructor%E_f%s_ptr%cells(i,j,k)
										constructor%v_s_f(dim,i,j,k) = constructor%v_s%s_ptr%cells(i,j,k)
										do spec = 1, species_number
											if (initial_riemann%get_success()) then
												constructor%Y_f(spec,dim,i,j,k) = Y_face_riemann(spec)
											else
												constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
											end if
										end do
										do dim1 = 1, dimensions
										if (dim1 /= dim) then
											constructor%v_f(dim1,dim,i,j,k) = constructor%v%v_ptr%pr(dim1)%cells(i,j,k)
										end if
										end do
									end if
								end if
							end if


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
		
		call constructor%state_eq%apply_state_equation_flow_variables_for_IC()

		call constructor%apply_boundary_conditions_main()
		
		! Same ambient-reference density convention as the FDS low-Mach solver.
		! It removes the hydrostatic reference force and retains density-driven
		! buoyancy in the compressible momentum and total-energy equations.
		constructor%rho_0 = constructor%rho%s_ptr%cells(cons_inner_loop(1,2), &
			cons_inner_loop(2,2),cons_inner_loop(3,2))
		if (constructor%rho_0 <= tiny(1.0_dp)) &
			error stop 'CABARET: non-positive buoyancy reference density'
		
        constructor%p_f     = constructor%p_f_new%s_ptr%cells
        constructor%rho_f   = constructor%rho_f_new%s_ptr%cells
		constructor%E_f_f	= constructor%E_f_f_new%s_ptr%cells
		constructor%v_s_f	= constructor%v_s_f_new%s_ptr%cells
       

        call manager%create_timer(cabaret_timer                 ,'CABARET solver time'              , 'sol_t')
        call manager%create_timer(cabaret_gas_dynamics_timer    ,'CABARET gas dynamics time'        , 'gd_t')
        call manager%create_timer(cabaret_eos_timer             ,'CABARET eos solver time'          , 'eos_t')
        call manager%create_timer(cabaret_chemistry_timer       ,'CABARET chemistry solver time'    , 'chem_t')
        call manager%create_timer(cabaret_diffusion_timer       ,'CABARET diffusion solver time'    , 'diff_t')
        call manager%create_timer(cabaret_heattransfer_timer    ,'CABARET heattransfer solver time' , 'ht_t')
        call manager%create_timer(cabaret_radiation_timer       ,'CABARET radiation solver time'    , 'rad_t')
        call manager%create_timer(cabaret_viscosity_timer       ,'CABARET viscosity solver time'    , 'visc_t')
	end function
	!> Advance the coupled CABARET solution by one full time step.
	subroutine solve_problem(this)
		class(cabaret_solver)	,intent(inout)	:: this

		call cabaret_timer%tic()

		! CABARET gas-dynamics stage: predictor -> face reconstruction -> corrector.
		! Physical source solvers are evaluated
		! once here.  Their full-step production fields are used as averaged
		! conservative source rates in both CABARET half-steps and in the
		! characteristic limiter shifts.
		call this%prepare_gas_dynamics_step()
		call this%evaluate_effective_physical_source_rates()

		call cabaret_gas_dynamics_timer%tic()
		call this%predict_conservative_half_step()
		call cabaret_gas_dynamics_timer%toc()

		call this%update_cell_thermodynamics()
		call cabaret_gas_dynamics_timer%tic()
		call this%reconstruct_acoustic_face_state()
		call this%reconstruct_contact_family_face_state()
		call this%finish_face_reconstruction()
		! Pressure and normal velocity must be equilibrated before density is
		! recomputed from the contact-family composition and temperature.
		call this%apply_material_contact_riemann_closure()
		call this%close_material_contact_density()
		call cabaret_gas_dynamics_timer%toc()

		call this%update_flow_thermodynamics()
		call cabaret_gas_dynamics_timer%tic()
		call this%correct_conservative_full_step()
		call this%finalize_gas_dynamics_step()
		call cabaret_gas_dynamics_timer%toc()

		! Dispersed-phase source rates are advanced and assembled together with
		! the other physical sources before the CABARET predictor.

		call this%update_cell_thermodynamics()
		call cabaret_gas_dynamics_timer%tic()
		call this%cache_flow_state_for_next_step()
		if (this%CFL_condition_flag) then
			call this%calculate_time_step()
		end if
		this%time = this%time + this%time_step
		call cabaret_gas_dynamics_timer%toc(new_iter=.true.)

		call cabaret_timer%toc(new_iter=.true.)
	end subroutine solve_problem


	!=======================================================================
	! CABARET gas-dynamics stages
	!=======================================================================

	!> Apply cell boundary conditions, exchange halos and cache the time-n state.

	subroutine prepare_gas_dynamics_step(this)
		class(cabaret_solver), intent(inout) :: this

		call this%apply_boundary_conditions_main()
		call this%exchange_conservative_state()
		call this%cache_conservative_state()
	end subroutine prepare_gas_dynamics_step

	!> Evaluate source solvers once and assemble conservative source rates.

	subroutine evaluate_effective_physical_source_rates(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, spec, phase
		integer, dimension(3,2) :: cons_inner_loop

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		! Use direct component association instead of assumed-shape dummy
		! arrays.  The effective-source arrays are allocated with conservative
		! ghost-cell bounds, e.g. 0:nx+1.  Passing them to a dummy declared as
		! dimension(:,:,:) would reset the lower bounds to 1 inside this
		! routine and shift all i,j,k indexing by one at boundaries.
		associate(  rho_src  => this%rho_src_old  , &
					mom_src  => this%mom_src_old  , &
					rhoE_src => this%rhoE_src_old , &
					rhoY_src => this%rhoY_src_old)

		rho_src  = 0.0_dp
		mom_src  = 0.0_dp
		rhoE_src = 0.0_dp
		rhoY_src = 0.0_dp


		! Evaluate selected physical source solvers once per full CABARET time step.
		! Existing production fields are interpreted as effective conservative source
		! rates, or as averaged rates if the solver internally integrates over dt.
		if (this%heat_trans_flag) then
			call cabaret_heattransfer_timer%tic()
			call this%heat_trans_solver%solve_heat_transfer(this%time_step)
			call cabaret_heattransfer_timer%toc(new_iter=.true.)
		end if

		if (this%radiation_flag) then
			call cabaret_radiation_timer%tic()
			call this%radiation_solver%solve_radiation(this%time_step)
			call cabaret_radiation_timer%toc(new_iter=.true.)
		end if
		if (this%diffusion_flag) then
			call cabaret_diffusion_timer%tic()
			call this%diff_solver%solve_diffusion(this%time_step)
			call cabaret_diffusion_timer%toc(new_iter=.true.)
		end if

		if (this%viscosity_flag) then
			call cabaret_viscosity_timer%tic()
			call this%viscosity_solver%solve_viscosity(this%time_step)
			call cabaret_viscosity_timer%toc(new_iter=.true.)
		end if

		if (this%reactive_flag) then
			call cabaret_chemistry_timer%tic()
			call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
			call cabaret_chemistry_timer%toc(new_iter=.true.)
		end if

		! Particle drag and heat transfer require current gas transport coefficients
		! even when viscous stress and Fourier conduction are disabled.
		if (this%additional_particles_phases_number > 0) then
			if (.not. this%heat_trans_flag) &
				call this%heat_trans_solver%update_thermal_conductivity()
			if (.not. this%viscosity_flag) &
				call this%viscosity_solver%update_dynamic_viscosity()
		end if

		! Advance exactly one runtime-selected dispersed-phase backend per phase.
		! The facade returns rates in the common gas-coupling contract.
		call this%advance_particle_phases()

		associate(  E_f_prod_chem  => this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat  => this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc  => this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff  => this%E_f_prod_diff%s_ptr	, &
					Y_prod_chem    => this%Y_prod_chem%v_ptr		, &
					Y_prod_diff    => this%Y_prod_diff%v_ptr		, &
					v_prod_visc    => this%v_prod_visc%v_ptr		, &
					bc             => this%boundary%bc_ptr)

		!$omp parallel default(shared) private(i,j,k,dim,spec,phase)
		!$omp do collapse(3) schedule(guided)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if (bc%bc_markers(i,j,k) == 0) then
				if (this%reactive_flag) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_chem%cells(i,j,k)
					do spec = 1, species_number
						rhoY_src(spec,i,j,k) = rhoY_src(spec,i,j,k) + Y_prod_chem%pr(spec)%cells(i,j,k)
					end do
				end if

				if (this%heat_trans_flag) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_heat%cells(i,j,k)
				end if

				if (this%radiation_flag) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + this%E_f_prod_rad%s_ptr%cells(i,j,k)
				end if
				if (this%diffusion_flag) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_diff%cells(i,j,k)
					do spec = 1, species_number
						rhoY_src(spec,i,j,k) = rhoY_src(spec,i,j,k) + Y_prod_diff%pr(spec)%cells(i,j,k)
					end do
				end if

				if (this%viscosity_flag) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_visc%cells(i,j,k)
					do dim = 1, dimensions
						mom_src(dim,i,j,k) = mom_src(dim,i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k)
					end do
				end if

				! Reduced-gravity buoyancy in the FDS velocity convention.  The
				! ambient hydrostatic contribution is represented by the reference
				! pressure; the remaining acceleration is (rho_0-rho)g/rho.
				do dim = 1, dimensions
					mom_src(dim,i,j,k) = mom_src(dim,i,j,k)+ &
						(this%rho_0-this%rho_old(i,j,k))*this%g(dim)
					rhoE_src(i,j,k) = rhoE_src(i,j,k)+ &
						(this%rho_0-this%rho_old(i,j,k))*this%v_old(dim,i,j,k)*this%g(dim)
				end do

				! Dispersed-phase source fields are rates.  The momentum field is
				! acceleration, so convert it to conservative momentum production.
				do phase = 1, this%additional_particles_phases_number
					rho_src(i,j,k) = rho_src(i,j,k) + &
						this%rho_prod_particles(phase)%s_ptr%cells(i,j,k)
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + &
						this%E_f_prod_particles(phase)%s_ptr%cells(i,j,k)
					do dim = 1, dimensions
						mom_src(dim,i,j,k) = mom_src(dim,i,j,k) + &
							this%rho_old(i,j,k) * &
							this%v_prod_particles(phase)%v_ptr%pr(dim)%cells(i,j,k)
					end do
					do spec = 1, species_number
						rhoY_src(spec,i,j,k) = rhoY_src(spec,i,j,k) + &
							this%Y_prod_particles(phase)%v_ptr%pr(spec)%cells(i,j,k)
					end do
				end do
			end if
		end do
		end do
		end do
		!$omp end do nowait
		!$omp end parallel
		end associate

		end associate

	end subroutine evaluate_effective_physical_source_rates


	!> Project conservative source rates onto the two acoustic characteristic families.


	subroutine physical_source_acoustic_shifts(this, dim, rho_state, p_state, c_state, E_state, vel_state, Y_state, &
			S_rho, S_mom, S_rhoE, S_rhoY, g_R, g_Q, g_contact, g_vel)
		class(cabaret_solver), intent(in) :: this
		integer, intent(in) :: dim
		real(dp), intent(in) :: rho_state, p_state, c_state, E_state
		real(dp), dimension(:), intent(in) :: vel_state, Y_state, S_mom, S_rhoY
		real(dp), intent(in) :: S_rho, S_rhoE
		real(dp), intent(out) :: g_R, g_Q, g_contact
		real(dp), dimension(:), intent(out) :: g_vel

		integer :: dim1
		real(dp) :: f_p, f_e, f_ud
		real(dp), dimension(size(Y_state)) :: g_Y_dummy

		! Exact thermally-perfect multicomponent pressure source projection.
		! The pressure rate is computed in the thermophysical module as
		!   dp/dt = p_rho*S_rho + p_eps*de/dt + sum_k p_Yk*dY_k/dt,
		! so chemistry and diffusion composition changes enter the acoustic source
		! shifts through the JANAF-based EOS, not through an effective-gamma closure.
		call this%thermo%thermo_ptr%pressure_source_rate_from_conservative_sources( &
			rho_state, p_state, E_state, vel_state, Y_state, S_rho, S_mom, S_rhoE, S_rhoY, &
			f_p, f_e, g_Y_dummy)

		f_ud = (S_mom(dim) - vel_state(dim)*S_rho) / rho_state
		g_R = f_ud + f_p/(rho_state*c_state)
		g_Q = f_ud - f_p/(rho_state*c_state)
		g_contact = f_p - c_state*c_state*S_rho

		do dim1 = 1, size(vel_state)
			g_vel(dim1) = (S_mom(dim1) - vel_state(dim1)*S_rho) / rho_state
		end do
	end subroutine physical_source_acoustic_shifts


	!> Project conservative source rates onto species and material temperature.
	!!
	!! Species rates follow dY_k/dt=(S_{rho Y_k}-Y_k S_rho)/rho.  The
	!! temperature rate is evaluated by a small thermodynamically consistent
	!! source probe through the JANAF mixture EOS.
	subroutine physical_source_contact_shifts(this, rho_state, p_state, E_state, vel_state, Y_state, &
			S_rho, S_mom, S_rhoE, S_rhoY, g_Y, g_temperature)
		class(cabaret_solver), intent(in) :: this
		real(dp), intent(in) :: rho_state, p_state, E_state
		real(dp), dimension(:), intent(in) :: vel_state, Y_state, S_mom, S_rhoY
		real(dp), intent(in) :: S_rho, S_rhoE
		real(dp), dimension(:), intent(out) :: g_Y
		real(dp), intent(out) :: g_temperature

		integer :: spec
		real(dp) :: pressure_rate, energy_rate, source_norm, probe_dt
		real(dp) :: rho_trial, p_trial, T_state, T_trial
		real(dp), dimension(size(Y_state)) :: Y_trial

		call this%thermo%thermo_ptr%pressure_source_rate_from_conservative_sources( &
			rho_state, p_state, E_state, vel_state, Y_state, S_rho, S_mom, S_rhoE, S_rhoY, &
			pressure_rate, energy_rate, g_Y)

		source_norm = abs(S_rho)/rho_state + abs(pressure_rate)/p_state
		do spec = 1, size(Y_state)
			source_norm = source_norm + abs(g_Y(spec))
		end do
		if (source_norm == 0.0_dp) then
			g_temperature = 0.0_dp
			return
		end if

		probe_dt = min(0.5_dp*this%time_step, 1.0e-8_dp/source_norm)
		rho_trial = rho_state + probe_dt*S_rho
		p_trial = p_state + probe_dt*pressure_rate
		do spec = 1, size(Y_state)
			Y_trial(spec) = Y_state(spec) + probe_dt*g_Y(spec)
		end do

		if (rho_trial <= 0.0_dp) error stop 'CABARET source probe: non-positive density'
		if (p_trial <= 0.0_dp) error stop 'CABARET source probe: non-positive pressure'
!		if (any(Y_trial < 0.0_dp)) error stop 'CABARET source probe: negative mass fraction'

		T_state = this%thermo%thermo_ptr%temperature_from_pressure_density_Y( &
			p_state, rho_state, Y_state)
		T_trial = this%thermo%thermo_ptr%temperature_from_pressure_density_Y( &
			p_trial, rho_trial, Y_trial)
		g_temperature = (T_trial - T_state)/probe_dt
	end subroutine physical_source_contact_shifts

	!> Close cell-centred thermodynamic variables with the mixture EOS.


	subroutine update_cell_thermodynamics(this)
		class(cabaret_solver), intent(inout) :: this

		call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation()
		call cabaret_eos_timer%toc(new_iter=.true.)
	end subroutine update_cell_thermodynamics

	!> Close face-centred thermodynamic variables with the mixture EOS.

	subroutine update_flow_thermodynamics(this)
		class(cabaret_solver), intent(inout) :: this

		call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation_flow_variables()
		call cabaret_eos_timer%toc(new_iter=.true.)
	end subroutine update_flow_thermodynamics

	!> Normalize and exchange reconstructed primitive face variables.

	subroutine finish_face_reconstruction(this)
		class(cabaret_solver), intent(inout) :: this

		call this%normalize_face_mass_fractions()
		call this%exchange_face_primitive_state()
	end subroutine finish_face_reconstruction
	!> Enforce pressure and normal-velocity equilibrium at material interfaces.
	!!
	!! The thermally-perfect Riemann solver supplies only the acoustic star state
	!! (p*,u*).  CABARET retains its contact-family mass fractions and material
	!! temperature; density is then recovered from rho=p*M(Y)/(R_u*T).
	subroutine apply_material_contact_riemann_closure(this)
		class(cabaret_solver), intent(inout) :: this

		type(riemann_solver) :: riemann
		integer :: dimensions, species_number
		integer :: i, j, k, dim, dim1, spec, il, jl, kl
		integer, dimension(3,2) :: face_inner, loop
		real(dp) :: molar_mass_left, molar_mass_right, molar_mass_ratio
		real(dp) :: rho_left, rho_right, p_left, p_right, u_left, u_right
		real(dp) :: gamma_left, gamma_right, c_left, c_right
		real(dp) :: p_floor, rho_floor, p_star, u_star, T_star, rho_star
		real(dp) :: p_scale, u_scale, T_scale, p_tol, u_tol, T_tol, face_molar_mass
		real(dp), dimension(:), allocatable :: Y_left, Y_right, Y_face
		logical :: invalid_face, material_face, equilibrium_contact

		dimensions = this%domain%get_domain_dimensions()
		species_number = this%chem%chem_ptr%species_number
		face_inner = this%domain%get_local_inner_faces_bounds()

		associate(rho => this%rho%s_ptr, p => this%p%s_ptr, T => this%T%s_ptr, &
			gamma => this%gamma%s_ptr, v => this%v%v_ptr, Y => this%Y%v_ptr, &
			rho_face => this%rho_f_new%s_ptr, p_face => this%p_f_new%s_ptr, &
			v_face => this%v_f_new%v_ptr, Y_face_field => this%Y_f_new%v_ptr, &
			bc => this%boundary%bc_ptr)
		!$omp parallel default(shared) private(riemann,i,j,k,dim,dim1,spec,il,jl,kl,loop, &
		!$omp& molar_mass_left,molar_mass_right,molar_mass_ratio,rho_left,rho_right,p_left,p_right, &
		!$omp& u_left,u_right,gamma_left,gamma_right,c_left,c_right,p_floor,rho_floor,p_star,u_star, &
		!$omp& T_star,rho_star,p_scale,u_scale,T_scale,p_tol,u_tol,T_tol,face_molar_mass, &
		!$omp& Y_left,Y_right,Y_face,invalid_face,material_face,equilibrium_contact)
		call riemann%clear(reset_counter=.true.)
		allocate(Y_left(species_number), Y_right(species_number), Y_face(species_number))

		do dim = 1, dimensions
			loop = face_inner
			do dim1 = 1, dimensions
				loop(dim1,2) = face_inner(dim1,2) - (1 - I_m(dim1,dim))
			end do

			!$omp do collapse(3) schedule(static)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				il = i - I_m(dim,1); jl = j - I_m(dim,2); kl = k - I_m(dim,3)
				if (bc%bc_markers(i,j,k) /= 0 .or. bc%bc_markers(il,jl,kl) /= 0) cycle

				molar_mass_left = 0.0_dp
				molar_mass_right = 0.0_dp
				do spec = 1, species_number
					Y_left(spec) = Y%pr(spec)%cells(il,jl,kl)
					Y_right(spec) = Y%pr(spec)%cells(i,j,k)
					Y_face(spec) = Y_face_field%pr(spec)%cells(dim,i,j,k)
					molar_mass_left = molar_mass_left + Y_left(spec)/this%thermo%thermo_ptr%molar_masses(spec)
					molar_mass_right = molar_mass_right + Y_right(spec)/this%thermo%thermo_ptr%molar_masses(spec)
				end do
				molar_mass_left = 1.0_dp/molar_mass_left
				molar_mass_right = 1.0_dp/molar_mass_right
				molar_mass_ratio = max(molar_mass_left,molar_mass_right)/min(molar_mass_left,molar_mass_right)

				invalid_face = p_face%cells(dim,i,j,k) <= 0.0_dp .or. rho_face%cells(dim,i,j,k) <= 0.0_dp
				material_face = molar_mass_ratio > material_contact_molar_mass_ratio
				if (.not. (invalid_face .or. material_face)) cycle

				rho_left = rho%cells(il,jl,kl); rho_right = rho%cells(i,j,k)
				p_left = p%cells(il,jl,kl); p_right = p%cells(i,j,k)
				u_left = v%pr(dim)%cells(il,jl,kl); u_right = v%pr(dim)%cells(i,j,k)
				gamma_left = gamma%cells(il,jl,kl); gamma_right = gamma%cells(i,j,k)
				p_floor = 1.0e-12_dp*max(1.0_dp,abs(p_left),abs(p_right))
				rho_floor = 1.0e-12_dp*max(1.0_dp,abs(rho_left),abs(rho_right))
				c_left = sqrt(max(gamma_left*p_left/max(rho_left,rho_floor),0.0_dp))
				c_right = sqrt(max(gamma_right*p_right/max(rho_right,rho_floor),0.0_dp))

				p_scale = max(1.0_dp,abs(p_left),abs(p_right))
				u_scale = max(1.0_dp,abs(u_left),abs(u_right),c_left,c_right)
				T_scale = max(1.0_dp,abs(T%cells(il,jl,kl)),abs(T%cells(i,j,k)))
				p_tol = max(contact_equilibrium_p_rel_tol*p_scale,1000.0_dp*epsilon(1.0_dp)*p_scale)
				u_tol = max(contact_equilibrium_u_rel_tol*u_scale,1000.0_dp*epsilon(1.0_dp)*u_scale)
				T_tol = max(contact_equilibrium_T_rel_tol*T_scale,1000.0_dp*epsilon(1.0_dp)*T_scale)
				equilibrium_contact = material_face .and. abs(p_left-p_right) <= p_tol .and. &
					abs(u_left-u_right) <= u_tol .and. abs(T%cells(il,jl,kl)-T%cells(i,j,k)) <= T_tol

				if (equilibrium_contact) then
					p_star = 0.5_dp*(p_left+p_right)
					u_star = 0.5_dp*(u_left+u_right)
					T_star = 0.5_dp*(T%cells(il,jl,kl)+T%cells(i,j,k))
				else
					call riemann%set_thermally_perfect_parameters(this%thermo%thermo_ptr, &
						rho_left,rho_right,p_left,p_right,u_left,u_right,Y_left,Y_right,p_floor,rho_floor)
					call riemann%solve()
					if (.not. riemann%get_success()) then
						if (invalid_face) error stop 'CABARET material-face Riemann solve failed'
						cycle
					end if
					p_star = riemann%get_pressure()
					u_star = riemann%get_velocity()
					T_star = this%material_temperature_f(dim,i,j,k)
				end if

				face_molar_mass = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_face)
				if (p_star <= 0.0_dp .or. T_star <= 0.0_dp) then
					if (invalid_face) error stop 'CABARET material-face closure produced invalid p or T'
					cycle
				end if
				rho_star = p_star*face_molar_mass/(r_gase_J*T_star)
				if (rho_star <= 0.0_dp) then
					if (invalid_face) error stop 'CABARET material-face closure produced invalid density'
					cycle
				end if

				p_face%cells(dim,i,j,k) = p_star
				rho_face%cells(dim,i,j,k) = rho_star
				v_face%pr(dim)%cells(dim,i,j,k) = u_star
				this%material_temperature_f(dim,i,j,k) = T_star
				this%material_molar_volume_f(dim,i,j,k) = face_molar_mass/rho_star
				this%material_sensible_energy_f(dim,i,j,k) = &
					this%thermo%thermo_ptr%sensible_energy_parameter_from_temperature_Y(T_star,Y_face)
			end do
			end do
			end do
			!$omp end do nowait
		end do

		deallocate(Y_left,Y_right,Y_face)
		!$omp end parallel
		end associate

		call this%exchange_face_primitive_state()
	end subroutine apply_material_contact_riemann_closure

	!> Conservative predictor from t^n to t^{n+1/2}.


	subroutine predict_conservative_half_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, spec
		integer, dimension(3,2) :: cons_inner_loop
		real(dp), dimension(3) :: cell_size
		character(len=20) :: coordinate_system
		real(dp) :: spec_summ, mean_higher, mean_lower, r, geom_coeff
		real(dp), dimension(:), allocatable :: rhoY_new

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		coordinate_system = this%domain%get_coordinate_system_name()
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		select case(coordinate_system)
		case ('cartesian')
			nu = 1
		case ('cylindrical')
			nu = 2
		case ('spherical')
			nu = 3
		case default
			nu = 1
		end select

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
					rho_src_old	=> this%rho_src_old, &
					mom_src_old	=> this%mom_src_old, &
					rhoE_src_old	=> this%rhoE_src_old, &
					rhoY_src_old	=> this%rhoY_src_old, &
					rho_rhs_old	=> this%rho_rhs_old, &
					mom_rhs_old	=> this%mom_rhs_old, &
					rhoE_rhs_old	=> this%rhoE_rhs_old, &
					rhoY_rhs_old	=> this%rhoY_rhs_old, &

					bc				=> this%boundary%bc_ptr		, &
					mesh			=> this%mesh%mesh_ptr)

        !$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r,geom_coeff,rhoY_new)
		allocate(rhoY_new(species_number))

		!$omp do collapse(3) schedule(guided)	 		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				
                r = mesh%mesh(1,i,j,k)
				geom_coeff = finite_volume_geometry_coefficient(nu, r, cell_size(1))
                
				rho%cells(i,j,k)	= 0.0_dp
				E_f%cells(i,j,k)	= 0.0_dp
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k)	=	0.0_dp 
				end do
		  
				do dim = 1,dimensions
					v%pr(dim)%cells(i,j,k)	=	0.0_dp 		
                end do				

				do dim = 1,dimensions
                    
                    mean_higher	= rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * v_f(dim,dim,i+i_m(dim,1), &
                    	& j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= rho_f(dim,i,j,k) * v_f(dim,dim,i,j,k)                   

					rho%cells(i,j,k)	=	rho%cells(i,j,k) - (mean_higher - mean_lower) /cell_size(1)
                    
                    if(dim == 1) then
                    	rho%cells(i,j,k) = rho%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                    end if
				end do
				
				rho_rhs_old(i,j,k) = rho%cells(i,j,k)
				rho%cells(i,j,k)		=	rho_old(i,j,k)	+  0.5_dp * this%time_step * &
					(rho_rhs_old(i,j,k) + rho_src_old(i,j,k))

				! Species are advanced in conservative partial-density form.
				! The stored primitive mass fractions are reconstructed only after the
				! rho*Y_k update, which keeps species transport tied to the mass flux.
				spec_summ = 0.0_dp
				do spec = 1,species_number
					rhoY_new(spec) = 0.0_dp
					do	dim = 1,dimensions
						mean_higher	=  rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * y_f(spec,dim,i+i_m(dim,1),j+i_m(dim,2), &
							& k+i_m(dim,3)) * v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f(dim,i,j,k) * y_f(spec,dim,i,j,k) * v_f(dim,dim,i,j,k)
						rhoY_new(spec)	=  rhoY_new(spec) - (mean_higher - mean_lower ) /cell_size(1)
						if(dim == 1) then
							rhoY_new(spec) = rhoY_new(spec) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
						end if
					end do
					rhoY_rhs_old(spec,i,j,k) = rhoY_new(spec)
					rhoY_new(spec) = rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * &
						(rhoY_rhs_old(spec,i,j,k) + rhoY_src_old(spec,i,j,k))
					spec_summ = spec_summ + max(rhoY_new(spec), 0.0_dp)
				end do

				if (spec_summ > 0.0_dp) then
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(rhoY_new(spec), 0.0_dp) / spec_summ
					end do
				else
					spec_summ = 0.0_dp
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(Y_old(spec,i,j,k), 0.0_dp)
						spec_summ = spec_summ + Y%pr(spec)%cells(i,j,k)
					end do
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k) / spec_summ
					end do
				end if
				do dim = 1,dimensions
					do dim1 = 1,dimensions
                        
						mean_higher	=  rho_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))*v_f(dim,dim1,i+I_m(dim1,1),j+I_m(dim1,2), &
							& k+I_m(dim1,3))*v_f(dim1,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))
						mean_lower	=  rho_f(dim1,i,j,k)*v_f(dim,dim1,i,j,k)*v_f(dim1,dim1,i,j,k)
                        
						v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower ) /cell_size(1)

                        if(dim1 == 1) then                        
                        	v%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                        end if                        

					end do

					v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k)	-	 ( p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p_f(dim, &
						& i,j,k)) /cell_size(1)
					
					mom_rhs_old(dim,i,j,k) = v%pr(dim)%cells(i,j,k)
					v%pr(dim)%cells(i,j,k)	=	rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dp*this%time_step * &
						(mom_rhs_old(dim,i,j,k) + mom_src_old(dim,i,j,k)) 
					v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k) /rho%cells(i,j,k)
	
				end do	
    
				do dim = 1,dimensions
                    
 					mean_higher	=  (rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*E_f_f(dim,i+I_m(dim,1),j+I_m(dim,2), &
 						& k+I_m(dim,3))	+	p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))*v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2), &
 						& k+I_m(dim,3))
					mean_lower	=  (rho_f(dim,i,j,k)*E_f_f(dim,i,j,k)																	+	p_f(dim,i,j,k))									*v_f(dim,dim,i,j,k)
                    
					E_f%cells(i,j,k)	= 	E_f%cells(i,j,k) -	(mean_higher - mean_lower ) /cell_size(1)
                    
                    if(dim == 1) then
                    	E_f%cells(i,j,k) = E_f%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                    end if
				end do	
				
				rhoE_rhs_old(i,j,k) = E_f%cells(i,j,k)
				E_f%cells(i,j,k)		=	rho_old(i,j,k) * E_f_old(i,j,k)  + 0.5_dp*this%time_step * &
					(rhoE_rhs_old(i,j,k) + rhoE_src_old(i,j,k))
				E_f%cells(i,j,k)		=	E_f%cells(i,j,k) / rho%cells(i,j,k) 
				
			end if
		end do
		end do
        end do
		!$omp end do nowait		
		deallocate(rhoY_new)
		!$omp end parallel				
        
                    end associate
	end subroutine predict_conservative_half_step


	!> Reconstruct pressure, density and velocity from acoustic quasi-invariants.


	subroutine reconstruct_acoustic_face_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, dim2, plus, sign, bound_number, spec
		integer, dimension(3,2) :: cons_inner_loop, cons_utter_loop, flow_inner_loop, loop
		real(dp), dimension(3) :: cell_size, characteristic_speed
		real(dp) :: mol_l_face, mol_r_face, mol_ratio_face
		logical :: material_contact_face, sound_point_lr, sound_point_rl
		character(len=20) :: coordinate_system
		real(dp), dimension(2) :: r_inv, q_inv, r_inv_new, q_inv_new
		real(dp), dimension(:,:), allocatable :: v_inv, v_inv_corrected, v_inv_new
		real(dp), dimension(:), allocatable :: v_inv_half, v_inv_old
		real(dp), dimension(:), allocatable :: vel_state, Y_state, S_mom, S_rhoY, g_vel
		real(dp) :: r_inv_half, q_inv_half, R_inv_old, Q_inv_old
		real(dp) :: G_half, G_half_old, G_half_lower, G_half_higher
		real(dp) :: g_inv, alpha_loc, diss_l, diss_r
		real(dp) :: geom_coeff, geom_source_R, geom_source_Q
		real(dp) :: phys_source_R, phys_source_Q, phys_source_contact
		real(dp) :: v_f_approx, v_s_f_approx

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		coordinate_system = this%domain%get_coordinate_system_name()
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		select case(coordinate_system)
		case ('cartesian')
			nu = 1
		case ('cylindrical')
			nu = 2
		case ('spherical')
			nu = 3
		case default
			nu = 1
		end select

		call this%exchange_conservative_state()

		call this%zero_new_flow_state()
        
        associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &

					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &

                    v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					p_f			=> this%p_f		, &
					Y_f			=> this%Y_f		, &
					
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					p_old		=> this%p_old	, &
					v_s_old		=> this%v_s_old	, &
					rho_src_old	=> this%rho_src_old, &
					mom_src_old	=> this%mom_src_old, &
					rhoE_src_old	=> this%rhoE_src_old, &
					rhoY_src_old	=> this%rhoY_src_old, &
					
					r_inv_corr      => this%r_inv_corr	, &
					q_inv_corr      => this%q_inv_corr	, &
					v_inv_corr      => this%v_inv_corr	, &            

					bc				=> this%boundary%bc_ptr, &
					mesh			=> this%mesh%mesh_ptr)
					
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,dim2,plus,loop,G_half,G_half_old,G_half_lower,G_half_higher, &
		!$omp& r_inv,R_inv_half,R_inv_old,q_inv,Q_inv_half,Q_inv_old,v_inv,v_inv_half,v_inv_old,r_inv_new,q_inv_new, &
		!$omp& v_inv_new,v_inv_corrected,g_inv,v_f_approx,v_s_f_approx,characteristic_speed,diss_l,diss_r,alpha_loc,sign, &
		!$omp& bound_number,spec,geom_coeff,geom_source_R,geom_source_Q,phys_source_R,phys_source_Q,phys_source_contact, &
		!$omp& vel_state,Y_state,S_mom,S_rhoY,g_vel)

		allocate(v_inv(dimensions,2), v_inv_new(dimensions,2))
		allocate(v_inv_half(dimensions), v_inv_old(dimensions))
		allocate(vel_state(dimensions), S_mom(dimensions), g_vel(dimensions))
		allocate(Y_state(species_number), S_rhoY(species_number))
     
		do dim = 1,dimensions
		
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
					r_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2), &
						& k+I_m(dim,3))
					R_inv_half	= v%pr(dim)%cells(i,j,k)								+ G_half*(p%cells(i,j,k))
					R_inv_old	= v_old(dim,i,j,k)										+ G_half_old*(p_old(i,j,k))
					
					q_inv(1) 	= v_f(dim,dim,i,j,k)									- G_half*p_f(dim,i,j,k)
					q_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2), &
						& k+I_m(dim,3))
					Q_inv_half	= v%pr(dim)%cells(i,j,k)								- G_half*(p%cells(i,j,k))
					Q_inv_old	= v_old(dim,i,j,k)										- G_half_old*(p_old(i,j,k))

					do dim1 = 1,dimensions
						if (dim1 == dim) then
							v_inv(dim,1)		= contact_quasi_invariant(p_f(dim,i,j,k), rho_f(dim,i,j,k), v_s%cells(i,j,k))
							v_inv(dim,2)		= contact_quasi_invariant(p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), &
											  rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)), v_s%cells(i,j,k)) 
							v_inv_half(dim)		= contact_quasi_invariant(p%cells(i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
							v_inv_old(dim)		= contact_quasi_invariant(p_old(i,j,k), rho_old(i,j,k), v_s_old(i,j,k))
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

					! Source-shifted maximum-principle interval.  Geometry is an in-step
					! balance source.  Physical sources use the same effective
					! conservative rates that have already entered the predictor and will
					! enter the corrector; no residual-based g_inv reconstruction is used.
					geom_coeff = 0.0_dp
					if (dim == 1) then
						geom_coeff = finite_volume_geometry_coefficient(nu, mesh%mesh(1,i,j,k), cell_size(1))
					end if
					geom_source_R = -geom_coeff * v_s%cells(i,j,k) * v%pr(1)%cells(i,j,k)
					geom_source_Q =  geom_coeff * v_s%cells(i,j,k) * v%pr(1)%cells(i,j,k)

					phys_source_R = 0.0_dp
					phys_source_Q = 0.0_dp
					phys_source_contact = 0.0_dp
					g_vel(:) = 0.0_dp
					do dim1 = 1, dimensions
						vel_state(dim1) = v%pr(dim1)%cells(i,j,k)
						S_mom(dim1) = mom_src_old(dim1,i,j,k)
					end do
					do spec = 1, species_number
						Y_state(spec) = Y%pr(spec)%cells(i,j,k)
						S_rhoY(spec) = rhoY_src_old(spec,i,j,k)
					end do
					call this%physical_source_acoustic_shifts(dim, rho%cells(i,j,k), p%cells(i,j,k), &
						v_s%cells(i,j,k), E_f%cells(i,j,k), vel_state, Y_state, &
						rho_src_old(i,j,k), S_mom, rhoE_src_old(i,j,k), S_rhoY, &
						phys_source_R, phys_source_Q, phys_source_contact, g_vel)

					g_inv = geom_source_R + phys_source_R
					r_inv_corr(1,dim,i,j,k) = limit_quasi_invariant(r_inv_new(1), r_inv(1), R_inv_half, r_inv(2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					r_inv_corr(2,dim,i,j,k) = limit_quasi_invariant(r_inv_new(2), r_inv(1), R_inv_half, r_inv(2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					
					g_inv = geom_source_Q + phys_source_Q
					q_inv_corr(1,dim,i,j,k) = limit_quasi_invariant(q_inv_new(1), q_inv(1), Q_inv_half, q_inv(2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					q_inv_corr(2,dim,i,j,k) = limit_quasi_invariant(q_inv_new(2), q_inv(1), Q_inv_half, q_inv(2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))

					do dim2 = 1,dimensions
						if (dim2 == dim) then
							g_inv = phys_source_contact
						else
							g_inv = g_vel(dim2)
						end if
						v_inv_corr(1,dim2,dim,i,j,k) = limit_quasi_invariant(v_inv_new(dim2,1), v_inv(dim2,1), v_inv_half(dim2), &
							& v_inv(dim2,2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
						v_inv_corr(2,dim2,dim,i,j,k) = limit_quasi_invariant(v_inv_new(dim2,2), v_inv(dim2,1), v_inv_half(dim2), &
							& v_inv(dim2,2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
		
                    end do
                    
                    !**************************** Boundary conditions *****************************
					if (( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,1) ) .and. &
						( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,2) )) then
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then

								v_f_approx		= 0.5_dp * (v%pr(dim)%cells(i,j,k)	+ v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2), &
									& k+sign*I_m(dim,3)))
								v_s_f_approx	= 0.5_dp * (v_s%cells(i,j,k)			+ v_s%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim, &
									& 3)))
		
								characteristic_speed(1) = v_f_approx + v_s_f_approx
								characteristic_speed(2) = v_f_approx - v_s_f_approx
								characteristic_speed(3) = v_f_approx		
									
								call this%apply_boundary_conditions_flow(dim, i,j,k, characteristic_speed, q_inv_corr(:,dim,i,j,k), &
									& r_inv_corr(:,dim,i,j,k), v_inv_corr(:,:,dim,i,j,k), G_half)
								
							end if
						end do
					end if
				end if
			end do
			end do
			end do
			!$omp end do
		
        end do

		deallocate(v_inv, v_inv_new, v_inv_half, v_inv_old)
		deallocate(vel_state, S_mom, g_vel, Y_state, S_rhoY)
		!$omp end parallel                    
        end associate
                    
        associate(	rho			=> this%rho%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
            
					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					
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
            
			!$omp parallel default(shared)  private(i,j,k,dim,dim1,loop,G_half_lower,G_half_higher, &
			!$omp& v_f_approx,v_s_f_approx,characteristic_speed,sign,bound_number,spec, &
			!$omp& mol_l_face,mol_r_face,mol_ratio_face,material_contact_face,sound_point_lr,sound_point_rl)

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

					G_half_lower	= 1.0_dp / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))*rho%cells(i-I_m(dim,1),j-I_m(dim,2), &
						& k-I_m(dim,3)))
					G_half_higher	= 1.0_dp / (v_s%cells(i,j,k)*rho%cells(i,j,k))
                        
					v_f_approx		= 0.5_dp*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v%pr(dim)%cells(i,j,k))
					v_s_f_approx	= 0.5_dp*(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))			+ v_s%cells(i,j,k))

						mol_l_face = 0.0_dp
						mol_r_face = 0.0_dp
						material_contact_face = .false.
						do spec = 1, species_number
							mol_l_face = mol_l_face + Y%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) / &
								this%thermo%thermo_ptr%molar_masses(spec)
							mol_r_face = mol_r_face + Y%pr(spec)%cells(i,j,k) / this%thermo%thermo_ptr%molar_masses(spec)
						end do
						if (mol_l_face > 0.0_dp .and. mol_r_face > 0.0_dp) then
							mol_l_face = 1.0_dp / mol_l_face
							mol_r_face = 1.0_dp / mol_r_face
							mol_ratio_face = max(mol_l_face,mol_r_face) / min(mol_l_face,mol_r_face)
							material_contact_face = (mol_ratio_face > material_contact_molar_mass_ratio)
						end if

						sound_point_lr = ((abs(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) < &
							abs(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))) .and. &
							(v%pr(dim)%cells(i,j,k) > v_s%cells(i,j,k)))
						sound_point_rl = ((v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) < &
							-v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) .and. &
							(abs(v%pr(dim)%cells(i,j,k)) < abs(v_s%cells(i,j,k))))
                    
						if	((material_contact_face .or. (.not. sound_point_lr)) .and. &
							(material_contact_face .or. (.not. sound_point_rl))) then 

						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx
					
						if (( characteristic_speed(1) >= 0.0_dp )	.and.&
							( characteristic_speed(2) < 0.0_dp )		.and.&
							( characteristic_speed(3) >= 0.0_dp )) then				

							p_f_new%cells(dim,i,j,k)			=	(r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	-   q_inv_corr(1,dim,i, &
								& j,k)) / (G_half_lower + G_half_higher)
							rho_f_new%cells(dim,i,j,k)			=	density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2, &
								& dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
											    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(G_half_higher * r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
										& 3)) + G_half_lower * q_inv_corr(1,dim,i,j,k))	/ (G_half_lower + G_half_higher)
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
						
							p_f_new%cells(dim,i,j,k)			=	(r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	-   q_inv_corr(1,dim,i, &
								& j,k)) / (G_half_lower + G_half_higher)
							rho_f_new%cells(dim,i,j,k)			=	density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1, &
								& dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
							
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(G_half_higher * r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
										& 3)) + G_half_lower * q_inv_corr(1,dim,i,j,k))	/ (G_half_lower + G_half_higher)
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
								end if
							end do
							continue
						end if

						if (( characteristic_speed(1) >= 0.0_dp ).and.&
							( characteristic_speed(2) >= 0.0_dp ).and.&
							( characteristic_speed(3) >= 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)	= 0.5_dp * (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) - q_inv_corr(2, &
								& dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) / G_half_lower
							rho_f_new%cells(dim,i,j,k)	= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim, &
								& dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
											    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
						
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	0.5_dp * (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
										& 3)) + q_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end if
							end do	
						end if
					
						if (( characteristic_speed(1) < 0.0_dp ).and.&
							( characteristic_speed(2) < 0.0_dp ).and.&
							( characteristic_speed(3) < 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)			= 0.5_dp * (r_inv_corr(1,dim,i,j,k) - q_inv_corr(1,dim,i,j,k)) / G_half_higher
							rho_f_new%cells(dim,i,j,k)			= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1, &
								& dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
						
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
     
						if ((.not. material_contact_face) .and. sound_point_lr) then 
						
                            
						v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.5_dp*(v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim, &
							& 1),j-I_m(dim,2),k-I_m(dim,3))/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) &
																*0.5_dp*(v_s%cells(i,j,k) + v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
		   
						p_f_new%cells(dim,i,j,k)			= (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim, &
							& 3)) - v_f_new%pr(dim)%cells(dim,i,j,k))/G_half_lower
     
                            
						if (characteristic_speed(3) >= 0.0_dp) then
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim, &
								& dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
										    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							
							do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim, &
								& dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
                                
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
                                end if
							end do
                        end if                            
						
                        continue
                            
					end if
     
						if ((.not. material_contact_face) .and. sound_point_rl) then 
					
						v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.5_dp*(v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) + v%pr(dim)%cells(i-I_m(dim, &
							& 1),j-I_m(dim,2),k-I_m(dim,3))/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))) &
																*0.5_dp*(v_s%cells(i,j,k) + v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
		   
						p_f_new%cells(dim,i,j,k)			= (v_f_new%pr(dim)%cells(dim,i,j,k) - q_inv_corr(1,dim,i,j,k))/G_half_higher
						
						if (characteristic_speed(3) >= 0.0_dp) then
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim, &
								& dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
										    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim, &
								& dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
                                
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
	end subroutine reconstruct_acoustic_face_state


	!> Reconstruct species, temperature, molar volume and sensible energy along u-characteristics.


	subroutine reconstruct_contact_family_face_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, spec, bound_number
		integer, dimension(3,2) :: cons_inner_loop, cons_utter_loop, loop
		real(dp), dimension(3) :: cell_size
		real(dp), dimension(:,:), allocatable :: Y_inv, Y_inv_corrected, Y_inv_new
		real(dp), dimension(:), allocatable :: Y_inv_half, Y_inv_old
		real(dp), dimension(:), allocatable :: vel_state, Y_state, S_mom, S_rhoY, g_Y
		real(dp) :: g_inv, alpha_loc, g_temperature, S_inv_lower_bound
		real(dp) :: v_f_approx_lower, v_f_approx_higher
		real(dp) :: diss_l, diss_r
		real(dp), dimension(2) :: S_inv, S_inv_new, S_inv_corrected
		real(dp), dimension(2) :: Psi_inv, Psi_inv_new, Psi_inv_corrected
		real(dp), dimension(2) :: Ksi_inv, Ksi_inv_new, Ksi_inv_corrected
		real(dp) :: S_inv_half, Psi_inv_half, Ksi_inv_half
		real(dp) :: T_face_lower, T_face_higher
		real(dp) :: m_face_lower, m_face_higher, m_cell, m_corr, molar_ratio_contact
		real(dp) :: y_pos, y_sum_corr, molar_denom_corr, progress

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		call this%exchange_face_pressure_density()
		call this%initialize_material_contact_faces()

		! Species are treated as material/contact-family invariants of the
		! inviscid multicomponent Euler subsystem, not as thermodynamically
		! passive scalars.  The final face choice is therefore made with the
		! reconstructed contact velocity v_f_new used by the gas-dynamic face
		! state, so rho_f, p_f and Y_f describe one coherent mixture state.

        associate(	v			=> this%v%v_ptr	            , &
					Y			=> this%Y%v_ptr	            , &
					T			=> this%T%s_ptr	            , &
					rho			=> this%rho%s_ptr	            , &
					p			=> this%p%s_ptr	            , &
					E_f			=> this%E_f%s_ptr	            , &
					v_s			=> this%v_s%s_ptr	            , &
                    Y_f			=> this%Y_f		            , &
					Y_old		=> this%Y_old	            , &
					rho_src_old	=> this%rho_src_old, &
					mom_src_old	=> this%mom_src_old, &
					rhoE_src_old	=> this%rhoE_src_old, &
					rhoY_src_old	=> this%rhoY_src_old, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					material_temperature_f => this%material_temperature_f, &
					material_molar_volume_f => this%material_molar_volume_f, &
					material_sensible_energy_f => this%material_sensible_energy_f, &
					bc			=> this%boundary%bc_ptr)
        
		!$omp parallel default(shared) &
		!$omp& private(i,j,k,dim,loop,spec,Y_inv,Y_inv_half,Y_inv_new,Y_inv_old, &
		!$omp& g_inv,alpha_loc,Y_inv_corrected,v_f_approx_lower,v_f_approx_higher, &
		!$omp& bound_number,diss_r,diss_l,S_inv,S_inv_new,S_inv_corrected,S_inv_half, &
		!$omp& Psi_inv,Psi_inv_new,Psi_inv_corrected,Psi_inv_half,Ksi_inv,Ksi_inv_new, &
		!$omp& Ksi_inv_corrected,Ksi_inv_half,T_face_lower,T_face_higher, &
		!$omp& m_face_lower,m_face_higher,m_cell,m_corr,molar_ratio_contact,y_pos, &
		!$omp& y_sum_corr,molar_denom_corr,progress,vel_state,Y_state,S_mom,S_rhoY,g_Y,g_temperature,S_inv_lower_bound)
		allocate(Y_inv(species_number,2), Y_inv_corrected(species_number,2), Y_inv_new(species_number,2))
		allocate(Y_inv_half(species_number), Y_inv_old(species_number))
		allocate(vel_state(dimensions), S_mom(dimensions))
		allocate(Y_state(species_number), S_rhoY(species_number), g_Y(species_number))
		
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
                    end do

					! Temperature, molar volume and sensible energy are reconstructed with the same
					! contact-family stencil as Y_k.  The production density closure uses
					! psi=M(Y)/rho, because rho=M/psi is the quantity that becomes
					! nonphysical when composition and thermodynamic contact layers drift.
					T_face_lower  = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(this%p_f(dim,i,j,k), &
						this%rho_f(dim,i,j,k), Y_f(:,dim,i,j,k))
					T_face_higher = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(this%p_f(dim,i+i_m(dim,1),j+i_m(dim, &
						& 2),k+i_m(dim,3)), &
						this%rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)), Y_f(:,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))
					S_inv(1) = T_face_lower
					S_inv(2) = T_face_higher
					S_inv_half = T%cells(i,j,k)

					m_face_lower  = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_f(:,dim,i,j,k))
					m_face_higher = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_f(:,dim,i+i_m(dim,1),j+i_m(dim, &
						& 2),k+i_m(dim,3)))
					m_cell        = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_inv_half(:))
					Psi_inv(1)   = m_face_lower/this%rho_f(dim,i,j,k)
					Psi_inv(2)   = m_face_higher/this%rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					Psi_inv_half = m_cell/rho%cells(i,j,k)

					Ksi_inv(1) = this%thermo%thermo_ptr%sensible_energy_parameter_from_temperature_Y( &
						T_face_lower, Y_f(:,dim,i,j,k))
					Ksi_inv(2) = this%thermo%thermo_ptr%sensible_energy_parameter_from_temperature_Y( &
						T_face_higher, Y_f(:,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))
					Ksi_inv_half = this%thermo%thermo_ptr%sensible_energy_parameter_from_temperature_Y( &
						T%cells(i,j,k), Y_inv_half(:))
					
					diss_l = 0.0_dp
					diss_r = 0.0_dp

					do spec = 1,species_number
						y_inv_new(spec,1)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_l)*y_inv(spec,2))/(1.0_dp+diss_l)
						y_inv_new(spec,2)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_r)*y_inv(spec,1))/(1.0_dp+diss_r)
					end do

					S_inv_new(1) = (2.0_dp*S_inv_half - (1.0_dp-diss_l)*S_inv(2))/(1.0_dp+diss_l)
					S_inv_new(2) = (2.0_dp*S_inv_half - (1.0_dp-diss_r)*S_inv(1))/(1.0_dp+diss_r)
					Psi_inv_new(1) = (2.0_dp*Psi_inv_half - (1.0_dp-diss_l)*Psi_inv(2))/(1.0_dp+diss_l)
					Psi_inv_new(2) = (2.0_dp*Psi_inv_half - (1.0_dp-diss_r)*Psi_inv(1))/(1.0_dp+diss_r)
					Ksi_inv_new(1) = (2.0_dp*Ksi_inv_half - (1.0_dp-diss_l)*Ksi_inv(2))/(1.0_dp+diss_l)
					Ksi_inv_new(2) = (2.0_dp*Ksi_inv_half - (1.0_dp-diss_r)*Ksi_inv(1))/(1.0_dp+diss_r)

					g_Y(:) = 0.0_dp
					g_temperature = 0.0_dp
					do spec = 1, species_number
						Y_state(spec) = Y%pr(spec)%cells(i,j,k)
						S_rhoY(spec) = rhoY_src_old(spec,i,j,k)
					end do
					do spec = 1, dimensions
						vel_state(spec) = v%pr(spec)%cells(i,j,k)
						S_mom(spec) = mom_src_old(spec,i,j,k)
					end do
					call this%physical_source_contact_shifts(rho%cells(i,j,k), p%cells(i,j,k), &
						E_f%cells(i,j,k), vel_state, Y_state, rho_src_old(i,j,k), S_mom, &
						rhoE_src_old(i,j,k), S_rhoY, g_Y, g_temperature)
					g_inv = g_temperature
					S_inv_lower_bound = 0.0_dp
					alpha_loc = 0.0_dp
					S_inv_corrected(1) = limit_quasi_invariant(S_inv_new(1), S_inv(1), S_inv_half, S_inv(2), &
						g_inv*this%time_step, alpha_loc, S_inv_lower_bound, huge(1.0_dp))
					S_inv_corrected(2) = limit_quasi_invariant(S_inv_new(2), S_inv(1), S_inv_half, S_inv(2), &
						g_inv*this%time_step, alpha_loc, S_inv_lower_bound, huge(1.0_dp))

					g_inv = 0.0_dp
					alpha_loc = 0.0_dp
					Psi_inv_corrected(1) = limit_quasi_invariant(Psi_inv_new(1), Psi_inv(1), Psi_inv_half, Psi_inv(2), &
						g_inv*this%time_step, alpha_loc, 0.0_dp, huge(1.0_dp))
					Psi_inv_corrected(2) = limit_quasi_invariant(Psi_inv_new(2), Psi_inv(1), Psi_inv_half, Psi_inv(2), &
						g_inv*this%time_step, alpha_loc, 0.0_dp, huge(1.0_dp))
					Ksi_inv_corrected(1) = limit_quasi_invariant(Ksi_inv_new(1), Ksi_inv(1), Ksi_inv_half, Ksi_inv(2), &
						g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					Ksi_inv_corrected(2) = limit_quasi_invariant(Ksi_inv_new(2), Ksi_inv(1), Ksi_inv_half, Ksi_inv(2), &
						g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					
					do spec = 1,species_number
						! Geometry changes rho and rho*Y_k by the same factor; physical
						! effective sources use dY_k/dt=(S_rhoY-Y_k*S_rho)/rho.
						g_inv = g_Y(spec)
						alpha_loc = 0.0_dp
						y_inv_corrected(spec,1) = limit_quasi_invariant(y_inv_new(spec,1), y_inv(spec,1), Y_inv_half(spec), y_inv(spec, &
							& 2), &
													  g_inv*this%time_step, alpha_loc, 0.0_dp, 1.0_dp)
						y_inv_corrected(spec,2) = limit_quasi_invariant(y_inv_new(spec,2), y_inv(spec,1), Y_inv_half(spec), y_inv(spec, &
							& 2), &
													  g_inv*this%time_step, alpha_loc, 0.0_dp, 1.0_dp)
					end do

					molar_ratio_contact = max(m_face_lower,m_face_higher)/min(m_face_lower,m_face_higher)
					if ((molar_ratio_contact > direct_psi_contact_molar_mass_ratio) .and. &
							(abs(m_face_higher - m_face_lower) > 10.0_dp*tiny(1.0_dp)*max(abs(m_face_higher),abs(m_face_lower)))) then
							do bound_number = 1, 2
								y_sum_corr = 0.0_dp
								molar_denom_corr = 0.0_dp
								do spec = 1, species_number
									y_pos = max(y_inv_corrected(spec,bound_number), 0.0_dp)
									y_sum_corr = y_sum_corr + y_pos
									molar_denom_corr = molar_denom_corr + y_pos/this%thermo%thermo_ptr%molar_masses(spec)
								end do
								if ((y_sum_corr > 0.0_dp) .and. (molar_denom_corr > 0.0_dp)) then
									m_corr = y_sum_corr/molar_denom_corr
									progress = (m_corr - m_face_lower)/(m_face_higher - m_face_lower)
									progress = min(1.0_dp, max(0.0_dp, progress))
									Psi_inv_corrected(bound_number) = Psi_inv(1) + progress*(Psi_inv(2) - Psi_inv(1))
								end if
							end do
						end if
					
#ifdef OMP					
					call omp_set_lock(lock(i,j,k))
					call omp_set_lock(lock(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))	
#endif

					!  Boundary conditions on Y are set up along with other flow variables, here values on boundaries should not be
					! updated
					do spec = 1,species_number
						if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,1) ) then
							bound_number	= bc%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
							if(bound_number == 0) then	
								! Consistency fix: upwind Y from the same reconstructed face velocity as the gas-dynamic state.
								v_f_approx_lower		= v_f_new%pr(dim)%cells(dim,i,j,k)
								if (v_f_approx_lower < 0.0_dp) then
                                	Y_f_new%pr(spec)%cells(dim,i,j,k) =  (y_inv_corrected(spec,1))
									if (spec == 1) then
										material_temperature_f(dim,i,j,k) = S_inv_corrected(1)
										material_molar_volume_f(dim,i,j,k) = Psi_inv_corrected(1)
										material_sensible_energy_f(dim,i,j,k) = Ksi_inv_corrected(1)
									end if
								end if	
							end if
						end if

						if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= cons_utter_loop(dim,2) ) then
							bound_number	= bc%bc_markers(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
							if(bound_number == 0) then
									! Consistency fix: use the reconstructed higher face velocity.
									v_f_approx_higher		= v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								if (v_f_approx_higher >= 0.0_dp) then
									Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (y_inv_corrected(spec,2)) 	
									if (spec == 1) then
										material_temperature_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = S_inv_corrected(2)
										material_molar_volume_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = Psi_inv_corrected(2)
										material_sensible_energy_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = Ksi_inv_corrected(2)
									end if
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
			
			! Face species are normalized once after the reconstruction loop,
			! in normalize_face_mass_fractions().

        end do
		deallocate(Y_inv, Y_inv_corrected, Y_inv_new, Y_inv_half, Y_inv_old)
		deallocate(vel_state, S_mom, Y_state, S_rhoY, g_Y)
		!$omp end parallel

        end associate
	end subroutine reconstruct_contact_family_face_state


	!> Conservative trapezoidal corrector from t^n to t^{n+1}.


	subroutine correct_conservative_full_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, spec
		integer, dimension(3,2) :: cons_inner_loop
		real(dp), dimension(3) :: cell_size
		character(len=20) :: coordinate_system
		real(dp) :: spec_summ, mean_higher, mean_lower, r, geom_coeff
		real(dp), dimension(:), allocatable :: rhoY_new

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		coordinate_system = this%domain%get_coordinate_system_name()
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		select case(coordinate_system)
		case ('cartesian')
			nu = 1
		case ('cylindrical')
			nu = 2
		case ('spherical')
			nu = 3
		case default
			nu = 1
		end select

		call this%exchange_face_thermodynamic_state()

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
					rho_src_old	=> this%rho_src_old, &
					mom_src_old	=> this%mom_src_old, &
					rhoE_src_old	=> this%rhoE_src_old, &
					rhoY_src_old	=> this%rhoY_src_old, &
					rho_rhs_old	=> this%rho_rhs_old, &
					mom_rhs_old	=> this%mom_rhs_old, &
					rhoE_rhs_old	=> this%rhoE_rhs_old, &
					rhoY_rhs_old	=> this%rhoY_rhs_old, &

					bc				=> this%boundary%bc_ptr     , &
					mesh			=> this%mesh%mesh_ptr)
            
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r,geom_coeff,rhoY_new)
		allocate(rhoY_new(species_number))
         
		!$omp do collapse(3) schedule(guided)			
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
  			
			if (bc%bc_markers(i,j,k) == 0) then	
                
                r = mesh%mesh(1,i,j,k)
				geom_coeff = finite_volume_geometry_coefficient(nu, r, cell_size(1))
                
				rho%cells(i,j,k)	= 0.0_dp
				E_f%cells(i,j,k)	= 0.0_dp
				
				do spec = 1,species_number
					Y%pr(spec)%cells(i,j,k)	=	0.0_dp 
				end do
		  
				do dim = 1,dimensions
					v%pr(dim)%cells(i,j,k)	=	0.0_dp 		
                end do	
                
				do dim = 1,dimensions
                    
                    mean_higher	= rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim, &
                    	& 3))  *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= rho_f_new%cells(dim,i,j,k)    *v_f_new%pr(dim)%cells(dim,i,j,k)                    

					rho%cells(i,j,k)	=	rho%cells(i,j,k) - (mean_higher - mean_lower) /cell_size(1)
                    
                    if(dim == 1) then
                    	rho%cells(i,j,k) = rho%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                    end if
				end do                

				rho%cells(i,j,k) = rho_old(i,j,k) + 0.5_dp*this%time_step * &
					(rho_rhs_old(i,j,k) + rho%cells(i,j,k)) + this%time_step*rho_src_old(i,j,k)
	   
				! Species are advanced as conservative partial densities rho*Y_k.
				! Mass fractions are recovered after positivity correction of rho*Y_k.
				spec_summ = 0.0_dp
				do spec = 1,species_number
					rhoY_new(spec) = 0.0_dp
					do	dim = 1,dimensions
						mean_higher	=  rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  * y_f_new%pr(spec)%cells(dim, &
							& i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f_new%cells(dim,i,j,k)    * y_f_new%pr(spec)%cells(dim,i,j,k)   *v_f_new%pr(dim)%cells(dim,i, &
							& j,k)
						rhoY_new(spec)	=  rhoY_new(spec) - (mean_higher - mean_lower ) /cell_size(1)
						if(dim == 1) then
							rhoY_new(spec) = rhoY_new(spec) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
						end if
					end do
					rhoY_new(spec) = rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * &
						(rhoY_rhs_old(spec,i,j,k) + rhoY_new(spec)) + this%time_step*rhoY_src_old(spec,i,j,k)
					spec_summ = spec_summ + max(rhoY_new(spec), 0.0_dp)
				end do

				if (spec_summ > 0.0_dp) then
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(rhoY_new(spec), 0.0_dp) / spec_summ
					end do
				else
					spec_summ = 0.0_dp
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = max(Y_old(spec,i,j,k), 0.0_dp)
						spec_summ = spec_summ + Y%pr(spec)%cells(i,j,k)
					end do
					do spec = 1,species_number
						Y%pr(spec)%cells(i,j,k) = Y%pr(spec)%cells(i,j,k) / spec_summ
					end do
				end if
				do dim = 1,dimensions
					do dim1 = 1,dimensions
						mean_higher	= rho_f_new%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))  *v_f_new%pr(dim)%cells(dim1, &
							& i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) *v_f_new%pr(dim1)%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2), &
							& k+i_m(dim1,3))
						mean_lower	= rho_f_new%cells(dim1,i,j,k)  *v_f_new%pr(dim)%cells(dim1,i,j,k) *v_f_new%pr(dim1)%cells(dim1,i,j,k)
													
						v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower)	/cell_size(1)
                        
                        
                        if(dim1 == 1) then                        
                        	v%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                        end if                        
					end do
	   
					mean_higher	= p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= p_f_new%cells(dim,i,j,k)
													
					v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(mean_higher - mean_lower)	/cell_size(1)
					
					v%pr(dim)%cells(i,j,k) = rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dp*this%time_step * &
						(mom_rhs_old(dim,i,j,k) + v%pr(dim)%cells(i,j,k)) + this%time_step*mom_src_old(dim,i,j,k)
					v%pr(dim)%cells(i,j,k)	= v%pr(dim)%cells(i,j,k) / rho%cells(i,j,k) 
				end do	
	   
				E_f%cells(i,j,k)		=	0.0_dp  
				do dim = 1,dimensions
					mean_higher	= (rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))*E_f_f_new%cells(dim,i+i_m(dim,1), &
						& j+i_m(dim,2),k+i_m(dim,3))	+	p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim, &
						& 3)))*v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
					mean_lower	= (rho_f_new%cells(dim,i,j,k)*E_f_f_new%cells(dim,i,j,k)	+	p_f_new%cells(dim,i,j, &
						& k))*v_f_new%pr(dim)%cells(dim,i,j,k)
				
					E_f%cells(i,j,k)        = 	E_f%cells(i,j,k)		-	(mean_higher - mean_lower)	/cell_size(1)
                    
                    if(dim == 1) then
                    	E_f%cells(i,j,k) = E_f%cells(i,j,k) - geom_coeff*0.5_dp*(mean_higher + mean_lower)
                    end if
				end do	
				E_f%cells(i,j,k) = rho_old(i,j,k) * E_f_old(i,j,k) + 0.5_dp*this%time_step * &
					(rhoE_rhs_old(i,j,k) + E_f%cells(i,j,k)) + this%time_step*rhoE_src_old(i,j,k)
				E_f%cells(i,j,k) = E_f%cells(i,j,k) /rho%cells(i,j,k)
	
			end if
        end do
		end do
		end do
        ! **************************************************************
		!$omp end do nowait		
		deallocate(rhoY_new)
		!$omp end parallel

        end associate
	end subroutine correct_conservative_full_step


	!> Apply final cell boundary conditions and exchange the corrected state.


	subroutine finalize_gas_dynamics_step(this)
		class(cabaret_solver), intent(inout) :: this

        associate(	v   => this%v%v_ptr)
			call this%mpi_support%exchange_conservative_vector_field(v)
        end associate

		call this%apply_boundary_conditions_main()
	end subroutine finalize_gas_dynamics_step
	!> Advance the configured dispersed-phase backend and refresh its source rates.
	subroutine advance_particle_phases(this)
		class(cabaret_solver), intent(inout) :: this
		integer :: phase

		do phase = 1, this%additional_particles_phases_number
			call this%particles_solver(phase)%advance(this%time_step, this%time)
		end do
	end subroutine advance_particle_phases

	!> Store reconstructed face variables as the old flow state for the next step.

	subroutine cache_flow_state_for_next_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, dim1, spec
		integer, dimension(3,2) :: flow_utter_loop

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		flow_utter_loop = this%domain%get_local_utter_faces_bounds()

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
            
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,spec)
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
	end subroutine cache_flow_state_for_next_step

	!=======================================================================
	! CABARET solver maintenance helpers
	!=======================================================================
	! The original implementation had several repeated OpenMP/MPI blocks inside
	! solve_problem().  The routines below collect the repeated state management
	! operations so that the numerical CABARET stages are easier to read and
	! modify.  They intentionally do not alter the mathematical algorithm except
	! for the explicitly marked consistency fixes in the face-species treatment.

	!> Cache the complete beginning-of-step cell state.

	subroutine cache_conservative_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, spec, dimensions, species_number

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number

		this%rho_old = this%rho%s_ptr%cells
		this%E_f_old = this%E_f%s_ptr%cells
		this%v_s_old = this%v_s%s_ptr%cells
		this%p_old   = this%p%s_ptr%cells

		do spec = 1, species_number
			this%Y_old(spec,:,:,:) = this%Y%v_ptr%pr(spec)%cells
		end do

		do dim = 1, dimensions
			this%v_old(dim,:,:,:) = this%v%v_ptr%pr(dim)%cells
		end do
	end subroutine cache_conservative_state

	!> Exchange all conservative/thermodynamic cell halo fields required by CABARET.

	subroutine exchange_conservative_state(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%rho%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%E_f%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%v_s%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%gamma%s_ptr)

		call this%mpi_support%exchange_conservative_vector_field(this%Y%v_ptr)
		call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)
	end subroutine exchange_conservative_state

	!> Exchange reconstructed face pressure and density.

	subroutine exchange_face_pressure_density(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
	end subroutine exchange_face_pressure_density

	!> Exchange reconstructed face pressure, density, velocity and composition.

	subroutine exchange_face_primitive_state(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_primitive_state

	!> Exchange all face fields needed by the conservative corrector.

	subroutine exchange_face_thermodynamic_state(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%e_i_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%E_f_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%v_s_f_new%s_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_thermodynamic_state

	!> Clear face reconstruction output fields before characteristic extrapolation.

	subroutine zero_new_flow_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, spec, dimensions, species_number

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number

		this%p_f_new%s_ptr%cells   = 0.0_dp
		this%rho_f_new%s_ptr%cells = 0.0_dp
		if (allocated(this%material_temperature_f)) this%material_temperature_f = 0.0_dp
		if (allocated(this%material_molar_volume_f)) this%material_molar_volume_f = 0.0_dp
		if (allocated(this%material_sensible_energy_f)) this%material_sensible_energy_f = 0.0_dp

		do dim = 1, dimensions
			this%v_f_new%v_ptr%pr(dim)%cells = 0.0_dp
		end do

		do spec = 1, species_number
			this%Y_f_new%v_ptr%pr(spec)%cells = 0.0_dp
		end do
	end subroutine zero_new_flow_state


	!> Initialize contact-family temperature, molar-volume and sensible-energy coordinates.
	!> Initialize contact-family temperature, molar volume and sensible energy.
	subroutine initialize_material_contact_faces(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, dim1, i, j, k, dimensions
		integer, dimension(3,2) :: face_inner, loop
		real(dp) :: T_face, molar_mass

		dimensions = this%domain%get_domain_dimensions()
		face_inner = this%domain%get_local_inner_faces_bounds()

		associate(p_face => this%p_f, rho_face => this%rho_f, Y_face => this%Y_f)
		!$omp parallel default(shared) private(dim,dim1,i,j,k,loop,T_face,molar_mass)
		do dim = 1, dimensions
			loop = face_inner
			do dim1 = 1, dimensions
				loop(dim1,2) = face_inner(dim1,2) - (1 - I_m(dim1,dim))
			end do

			!$omp do collapse(3) schedule(static)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				T_face = this%thermo%thermo_ptr%temperature_from_pressure_density_Y( &
					p_face(dim,i,j,k), rho_face(dim,i,j,k), Y_face(:,dim,i,j,k))
				molar_mass = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions( &
					Y_face(:,dim,i,j,k))
				this%material_temperature_f(dim,i,j,k) = T_face
				this%material_molar_volume_f(dim,i,j,k) = molar_mass/rho_face(dim,i,j,k)
				this%material_sensible_energy_f(dim,i,j,k) = &
					this%thermo%thermo_ptr%sensible_energy_parameter_from_temperature_Y( &
						T_face, Y_face(:,dim,i,j,k))
			end do
			end do
			end do
			!$omp end do nowait
		end do
		!$omp end parallel
		end associate
	end subroutine initialize_material_contact_faces

	!> Close material-face density and enforce a local thermodynamic invariant domain.
	!!
	!! Baseline density is obtained from reconstructed (p,Y,T).  A selectively
	!! reconstructed molar-volume state and a pressure-equilibrium-preserving
	!! sensible-energy state may replace the baseline only inside the mixed part
	!! of a material layer.  The final density is bounded by adjacent cell states.
	subroutine close_material_contact_density(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: dim, dim1, i, j, k, il, jl, kl, spec
		integer, dimension(3,2) :: face_inner, loop
		real(dp) :: M_left, M_right, M_min, M_max, M_span, M_face, M_margin, M_ratio
		real(dp) :: T_left, T_right, T_guess, T_contact, T_reference
		real(dp) :: psi_used, psi_candidate, T_candidate, T_candidate_limit
		real(dp) :: psi_left, psi_right, psi_local_low, psi_local_high
		real(dp) :: ksi_candidate, T_base, T_pep, psi_pep
		real(dp) :: T_low, T_high, psi_T_low, psi_T_high
		real(dp) :: rho_min, rho_max, rho_span, rho_tol, rho_low, rho_high
		real(dp) :: psi_rho_low, psi_rho_high, psi_adm_low, psi_adm_high
		real(dp) :: dpsi, theta_low, theta_high, theta, psi_test, rho_face_value
		real(dp) :: molar_denom_left, molar_denom_right
		real(dp), dimension(:), allocatable :: Y_face
		logical :: mixed_face, direct_psi_valid

		dimensions = this%domain%get_domain_dimensions()
		species_number = this%chem%chem_ptr%species_number
		face_inner = this%domain%get_local_inner_faces_bounds()

		associate(rho_face => this%rho_f_new%s_ptr, p_face => this%p_f_new%s_ptr, &
			Y_face_field => this%Y_f_new%v_ptr, Y => this%Y%v_ptr, rho => this%rho%s_ptr, &
			T => this%T%s_ptr, bc => this%boundary%bc_ptr)
		!$omp parallel default(shared) private(dim,dim1,i,j,k,il,jl,kl,spec,loop,Y_face, &
		!$omp& M_left,M_right,M_min,M_max,M_span,M_face,M_margin,M_ratio,T_left,T_right, &
		!$omp& T_guess,T_contact,T_reference,psi_used,psi_candidate,T_candidate,T_candidate_limit, &
		!$omp& psi_left,psi_right,psi_local_low,psi_local_high,ksi_candidate,T_base,T_pep,psi_pep,T_low,T_high, &
		!$omp& psi_T_low,psi_T_high,rho_min,rho_max,rho_span,rho_tol,rho_low, &
		!$omp& rho_high,psi_rho_low,psi_rho_high,psi_adm_low,psi_adm_high,dpsi,theta_low,theta_high, &
		!$omp& theta,psi_test,rho_face_value,molar_denom_left,molar_denom_right,mixed_face,direct_psi_valid)
		allocate(Y_face(species_number))

		do dim = 1, dimensions
			loop = face_inner
			do dim1 = 1, dimensions
				loop(dim1,2) = face_inner(dim1,2) - (1 - I_m(dim1,dim))
			end do

			!$omp do collapse(3) schedule(guided)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				il = i-I_m(dim,1); jl = j-I_m(dim,2); kl = k-I_m(dim,3)
				if (bc%bc_markers(i,j,k) /= 0 .or. bc%bc_markers(il,jl,kl) /= 0) cycle

				molar_denom_left = 0.0_dp
				molar_denom_right = 0.0_dp
				do spec = 1, species_number
					molar_denom_left = molar_denom_left + &
						Y%pr(spec)%cells(il,jl,kl)/this%thermo%thermo_ptr%molar_masses(spec)
					molar_denom_right = molar_denom_right + &
						Y%pr(spec)%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec)
					Y_face(spec) = Y_face_field%pr(spec)%cells(dim,i,j,k)
				end do
				M_left = 1.0_dp/molar_denom_left
				M_right = 1.0_dp/molar_denom_right
				M_min = min(M_left,M_right); M_max = max(M_left,M_right)
				M_ratio = M_max/M_min
				if (M_ratio <= contact_reconstruction_molar_mass_ratio) cycle

				T_contact = this%material_temperature_f(dim,i,j,k)
				if (T_contact <= 0.0_dp) error stop 'CABARET contact closure: non-positive material temperature'
				M_face = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_face)
				psi_used = r_gase_J*T_contact/p_face%cells(dim,i,j,k)

				T_left = T%cells(il,jl,kl); T_right = T%cells(i,j,k)
				T_guess = this%thermo%thermo_ptr%temperature_from_pressure_density_Y( &
					p_face%cells(dim,i,j,k),rho_face%cells(dim,i,j,k),Y_face)
				T_reference = max(T_guess,T_left,T_right)
				M_span = M_max-M_min

				! Select direct CABARET molar volume only in the genuinely mixed layer.
				M_margin = direct_psi_face_mixture_rel_tol*M_span
				mixed_face = M_ratio > direct_psi_contact_molar_mass_ratio .and. &
					M_span > 10.0_dp*tiny(1.0_dp)*M_max .and. &
					M_face > M_min+M_margin .and. M_face < M_max-M_margin
				psi_candidate = this%material_molar_volume_f(dim,i,j,k)
				T_candidate = p_face%cells(dim,i,j,k)*psi_candidate/r_gase_J
				T_candidate_limit = max(T_reference,T_contact)*(1.0_dp+direct_psi_temperature_rel_tol) + &
					direct_psi_temperature_abs_tol
				direct_psi_valid = mixed_face .and. psi_candidate > 0.0_dp .and. &
					T_candidate <= T_candidate_limit
				if (direct_psi_valid) psi_used = psi_candidate

				! Parameter-free sensible-energy PEP projection onto the local
				! temperature-density-molar-volume invariant domain.
				mixed_face = M_span > 10.0_dp*tiny(1.0_dp)*M_max .and. &
					M_face > M_min .and. M_face < M_max
				if (mixed_face) then
					ksi_candidate = this%material_sensible_energy_f(dim,i,j,k)
					T_base = p_face%cells(dim,i,j,k)*psi_used/r_gase_J
					T_pep = this%thermo%thermo_ptr%temperature_from_sensible_energy_parameter_Y( &
						ksi_candidate,Y_face,T_base)
					if (T_pep > 0.0_dp) then
						psi_pep = r_gase_J*T_pep/p_face%cells(dim,i,j,k)
					else
						psi_pep = -1.0_dp
					end if

					T_low = min(T_left,T_right,T_base)
					T_high = max(T_left,T_right,T_base)
					psi_T_low = r_gase_J*T_low/p_face%cells(dim,i,j,k)
					psi_T_high = r_gase_J*T_high/p_face%cells(dim,i,j,k)

					rho_min = min(rho%cells(il,jl,kl),rho%cells(i,j,k))
					rho_max = max(rho%cells(il,jl,kl),rho%cells(i,j,k))
					rho_span = rho_max-rho_min
					rho_tol = face_density_guard_abs_tol+face_density_guard_rel_tol*rho_span
					rho_low = rho_min-rho_tol; rho_high = rho_max+rho_tol
					if (rho_low > 0.0_dp) then
						psi_rho_low = M_face/rho_high
						psi_rho_high = M_face/rho_low
					else
						psi_rho_low = huge(1.0_dp); psi_rho_high = -huge(1.0_dp)
					end if

					psi_left = M_left/rho%cells(il,jl,kl)
					psi_right = M_right/rho%cells(i,j,k)
					psi_local_low = min(psi_left,psi_right)
					psi_local_high = max(psi_left,psi_right)
					psi_adm_low = max(psi_T_low,psi_rho_low,psi_local_low)
					psi_adm_high = min(psi_T_high,psi_rho_high,psi_local_high)

					dpsi = psi_pep-psi_used
					theta_low = 0.0_dp; theta_high = 1.0_dp
					if (psi_pep > 0.0_dp .and. psi_adm_low <= psi_adm_high .and. &
						abs(dpsi) > 10.0_dp*tiny(1.0_dp)*max(abs(psi_used),abs(psi_pep))) then
						if (dpsi > 0.0_dp) then
							theta_low = max(theta_low,(psi_adm_low-psi_used)/dpsi)
							theta_high = min(theta_high,(psi_adm_high-psi_used)/dpsi)
						else
							theta_low = max(theta_low,(psi_adm_high-psi_used)/dpsi)
							theta_high = min(theta_high,(psi_adm_low-psi_used)/dpsi)
						end if
						theta_low = max(0.0_dp,theta_low); theta_high = min(1.0_dp,theta_high)
						if (theta_high >= theta_low .and. theta_high > 0.0_dp) then
							theta = theta_high
							psi_test = psi_used+theta*dpsi
							if (psi_test >= psi_adm_low .and. psi_test <= psi_adm_high) psi_used = psi_test
						end if
					end if
				end if

				! Final density invariant-domain projection.
				if (M_ratio > face_density_guard_molar_mass_ratio) then
					rho_min = min(rho%cells(il,jl,kl),rho%cells(i,j,k))
					rho_max = max(rho%cells(il,jl,kl),rho%cells(i,j,k))
					rho_span = rho_max-rho_min
					rho_tol = face_density_guard_abs_tol+face_density_guard_rel_tol*rho_span
					rho_low = rho_min-rho_tol; rho_high = rho_max+rho_tol
					if (rho_low > 0.0_dp) then
						rho_face_value = M_face/psi_used
						if (rho_face_value < rho_low) psi_used = M_face/rho_low
						if (rho_face_value > rho_high) psi_used = M_face/rho_high
					end if
				end if
				rho_face%cells(dim,i,j,k) = M_face/psi_used
			end do
			end do
			end do
			!$omp end do nowait
		end do

		deallocate(Y_face)
		!$omp end parallel
		end associate

		call this%exchange_face_pressure_density()
	end subroutine close_material_contact_density

	!> Enforce non-negative face mass fractions with unit sum.

	subroutine normalize_face_mass_fractions(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, dim1, spec, i, j, k
		integer :: dimensions, species_number
		integer, dimension(3,2) :: flow_inner_loop, loop
		real(dp) :: spec_summ

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()

		do dim = 1, dimensions
			loop = flow_inner_loop
			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				if ((this%boundary%bc_ptr%bc_markers(i,j,k) == 0) .or. &
					(this%boundary%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) == 0)) then

					! First try the freshly reconstructed face mass fractions.
					spec_summ = 0.0_dp
					do spec = 1, species_number
						this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k), 0.0_dp)
						spec_summ = spec_summ + this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
					end do

					! If a boundary or a degenerate characteristic decision left the face unset,
					! fall back to the old face state.  This avoids a zero molecular weight in
					! apply_state_equation_flow_variables().
					if (spec_summ <= tiny(1.0_dp)) then
						spec_summ = 0.0_dp
						do spec = 1, species_number
							this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%Y_f(spec,dim,i,j,k), 0.0_dp)
							spec_summ = spec_summ + this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
						end do
					end if

					! Last-resort fallback: use the right cell state if available, otherwise
					! the left cell state.  This branch should only be hit for initially bad
					! or boundary-adjacent faces.
					if (spec_summ <= tiny(1.0_dp)) then
						spec_summ = 0.0_dp
						if (this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then
							do spec = 1, species_number
								this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%Y%v_ptr%pr(spec)%cells(i,j,k), 0.0_dp)
								spec_summ = spec_summ + this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
							end do
						else
							do spec = 1, species_number
								this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2), &
									& k-I_m(dim,3)), 0.0_dp)
								spec_summ = spec_summ + this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
							end do
						end if
					end if

					if (spec_summ > tiny(1.0_dp)) then
						do spec = 1, species_number
							this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) / spec_summ
						end do
					end if
				end if
			end do
			end do
			end do
		end do
	end subroutine normalize_face_mass_fractions


	subroutine apply_boundary_conditions_flow(this, dim,i,j,k, characteristic_speed, q_inv_corrected, r_inv_corrected, &
		& v_inv_corrected, G_half)

		class(cabaret_solver)		,intent(inout)		:: this
		integer						,intent(in)			:: i, j, k, dim
		real(dp)	,dimension(3)	,intent(in)			:: characteristic_speed
		real(dp)	,dimension(2)	,intent(in)			:: q_inv_corrected, r_inv_corrected
		real(dp)	,dimension(:,:)	,intent(in)			:: v_inv_corrected
		real(dp)					,intent(in)			:: G_half
		
		real(dp)				:: r_inf, q_inf, s_inf, G_half_inf

		real(dp)				:: spec_summ
		real(dp)				:: n_out, denom, r_face, q_face, s_face, p_face, u_face
		real(dp)				:: G_r, G_q, rho_ref, c_ref
		real(dp)				:: farfield_density, farfield_pressure, farfield_velocity
		real(dp)				:: p_inf, rho_inf, c_inf, u_inf

		real(dp)	,dimension(3)	:: cell_size
		character(len=20)		:: boundary_type_name
		integer					:: dimensions, species_number
		character(len=20)		:: coordinate_system
		integer					:: sign, bound_number, side
		integer					:: face_i, face_j, face_k, ghost_i, ghost_j, ghost_k
		integer 				:: plus, dim1, spec	
		logical				:: r_out, q_out, contact_out

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
						! Slip reflective face state for solid walls and symmetry planes.
						! This branch is intentionally independent of the local Mach/sign pattern:
						! the impermeability condition sets the normal face velocity to zero.
						if (sign == -1) then
							side   = 1
							face_i = i
							face_j = j
							face_k = k
							p_face = -q_inv_corrected(side)/G_half
						else
							side   = 2
							face_i = i + I_m(dim,1)
							face_j = j + I_m(dim,2)
							face_k = k + I_m(dim,3)
							p_face = r_inv_corrected(side)/G_half
						end if

						p_f_new%cells(dim,face_i,face_j,face_k) = p_face
						rho_f_new%cells(dim,face_i,face_j,face_k) = &
							density_from_contact_quasi_invariant(p_face, v_inv_corrected(side,dim), &
								rho%cells(i,j,k), v_s%cells(i,j,k))
						v_f_new%pr(dim)%cells(dim,face_i,face_j,face_k) = 0.0_dp
						do dim1 = 1, dimensions
							if (dim1 /= dim) then
								v_f_new%pr(dim1)%cells(dim,face_i,face_j,face_k) = v_inv_corrected(side,dim1)
							end if
						end do

						spec_summ = 0.0_dp
						do spec = 1, species_number
							Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = max(Y%pr(spec)%cells(i,j,k), 0.0_dp)
							spec_summ = spec_summ + Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k)
						end do
						if (spec_summ > 0.0_dp) then
							do spec = 1, species_number
								Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = &
									Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) / spec_summ
							end do
						end if

					case ('outlet','inlet')
						! Characteristic open-boundary treatment written in outward-normal form.
						! For each characteristic family, sign*lambda > 0 means the wave leaves
						! the computational domain and the interior CABARET extrapolated value is
						! used. sign*lambda < 0 means the wave enters the domain and the
						! prescribed/far-field state is used.  This state is assembled here from
						! the boundary object, not from the transport/source ghost cell: outlet
						! ghost cells may be zero-gradient for viscosity/heat/diffusion.

						if (sign == -1) then
							side   = 1
							face_i = i
							face_j = j
							face_k = k
						else
							side   = 2
							face_i = i + I_m(dim,1)
							face_j = j + I_m(dim,2)
							face_k = k + I_m(dim,3)
						end if

						ghost_i = i + sign*I_m(dim,1)
						ghost_j = j + sign*I_m(dim,2)
						ghost_k = k + sign*I_m(dim,3)
						n_out  = real(sign, dp)

						farfield_pressure = bc%boundary_types(bound_number)%get_farfield_pressure()
						farfield_density  = bc%boundary_types(bound_number)%get_farfield_density()
						farfield_velocity = bc%boundary_types(bound_number)%get_farfield_velocity()

						p_inf   = p%cells(ghost_i,ghost_j,ghost_k)
						rho_inf = rho%cells(ghost_i,ghost_j,ghost_k)
						c_inf   = v_s%cells(ghost_i,ghost_j,ghost_k)
						u_inf   = farfield_velocity
						if (farfield_pressure > 0.0_dp) p_inf = farfield_pressure
						if (farfield_density  > 0.0_dp) rho_inf = farfield_density
						if (rho_inf <= 0.0_dp) error stop 'CABARET open boundary: non-positive density'
						if (p_inf   <= 0.0_dp) error stop 'CABARET open boundary: non-positive pressure'
						if (c_inf <= 0.0_dp) c_inf = v_s%cells(i,j,k)
						if (c_inf <= 0.0_dp) error stop 'CABARET open boundary: non-positive sound speed'

						G_half_inf = 1.0_dp / (rho_inf*c_inf)
						r_inf = u_inf + G_half_inf * p_inf
						q_inf = u_inf - G_half_inf * p_inf
						s_inf = contact_quasi_invariant(p_inf, rho_inf, c_inf)

						r_out       = (n_out*characteristic_speed(1) >= 0.0_dp)
						q_out       = (n_out*characteristic_speed(2) >= 0.0_dp)
						contact_out = (n_out*characteristic_speed(3) >= 0.0_dp)

						if (r_out) then
							r_face = r_inv_corrected(side)
							G_r    = G_half
						else
							r_face = r_inf
							G_r    = G_half_inf
						end if

						if (q_out) then
							q_face = q_inv_corrected(side)
							G_q    = G_half
						else
							q_face = q_inf
							G_q    = G_half_inf
						end if

						denom = G_r + G_q
						p_face = (r_face - q_face)/denom
						u_face = (G_q*r_face + G_r*q_face)/denom

						if (contact_out) then
							s_face  = v_inv_corrected(side,dim)
							rho_ref = rho%cells(i,j,k)
							c_ref   = v_s%cells(i,j,k)
						else
							s_face  = s_inf
							rho_ref = rho_inf
							c_ref   = c_inf
						end if

						p_f_new%cells(dim,face_i,face_j,face_k)   = p_face
						rho_f_new%cells(dim,face_i,face_j,face_k) = &
							density_from_contact_quasi_invariant(p_face, s_face, rho_ref, c_ref)

						do dim1 = 1, dimensions
							if (dim1 == dim) then
								v_f_new%pr(dim1)%cells(dim,face_i,face_j,face_k) = u_face
							else
								if (contact_out) then
									v_f_new%pr(dim1)%cells(dim,face_i,face_j,face_k) = v_inv_corrected(side,dim1)
								else
									v_f_new%pr(dim1)%cells(dim,face_i,face_j,face_k) = v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k)
								end if
							end if
						end do

						spec_summ = 0.0_dp
						do spec = 1, species_number
							if (contact_out) then
								Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = Y%pr(spec)%cells(i,j,k)
							else
								Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = Y%pr(spec)%cells(ghost_i,ghost_j,ghost_k)
							end if
							Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = &
								max(Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k), 0.0_dp)
							spec_summ = spec_summ + Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k)
						end do
						if (spec_summ > 0.0_dp) then
							do spec = 1, species_number
								Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = &
									Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) / spec_summ
							end do
						else
							do spec = 1, species_number
								Y_f_new%pr(spec)%cells(dim,face_i,face_j,face_k) = Y%pr(spec)%cells(i,j,k)
							end do
						end if

					case default
						write(*,*) 'CABARET apply_boundary_conditions_flow: unsupported boundary type ', trim(boundary_type_name)
						stop

				end select
			end if
		end do
		end associate

	end subroutine

	!> Fill conservative ghost cells for wall, symmetry and open boundaries.

	subroutine apply_boundary_conditions_main(this)

		class(cabaret_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer	,dimension(3,2)	:: cons_utter_loop, cons_inner_loop
		character(len=20)		:: boundary_type_name
		real(dp)				:: farfield_density, farfield_pressure, farfield_temperature, wall_temperature, farfield_velocity
		real(dp)				:: spec_summ, mol_mix_from_farfield, normal_velocity

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
		integer :: ghost_i, ghost_j, ghost_k

		dimensions			= this%domain%get_domain_dimensions()
		cons_utter_loop		= this%domain%get_local_utter_cells_bounds()	
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()	
        

		associate(  T				=> this%T%s_ptr					, &
					mix_mol_mass	=> this%mix_mol_mass%s_ptr		, &
					p				=> this%p%s_ptr					, &
					rho				=> this%rho%s_ptr				, &
					v				=> this%v%v_ptr					, &
					v_s				=> this%v_s%s_ptr				, &
					Y				=> this%Y%v_ptr					, &
					bc				=> this%boundary%bc_ptr)

		!$omp parallel default(shared)  private(i,j,k,plus,dim,dim1,sign,bound_number, &
		!$omp& ghost_i,ghost_j,ghost_k,farfield_pressure,farfield_density,farfield_temperature, &
		!$omp& farfield_velocity,mol_mix_from_farfield,spec_summ,normal_velocity, &
		!$omp& wall_temperature,boundary_type_name,specie_number)
		!$omp do collapse(3) schedule(static)

			do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
			do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
			do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,dimensions
						do plus = 1,2
							sign			= (-1)**plus
							ghost_i = i + sign*I_m(dim,1)
							ghost_j = j + sign*I_m(dim,2)
							ghost_k = k + sign*I_m(dim,3)

							bound_number	= bc%bc_markers(ghost_i,ghost_j,ghost_k)
							if( bound_number /= 0 ) then

								boundary_type_name = bc%boundary_types(bound_number)%get_type_name()
								select case(boundary_type_name)
									case('wall')

										p%cells(ghost_i,ghost_j,ghost_k)					= p%cells(i,j,k)
										rho%cells(ghost_i,ghost_j,ghost_k)				= rho%cells(i,j,k)
										T%cells(ghost_i,ghost_j,ghost_k)					= T%cells(i,j,k)
										mix_mol_mass%cells(ghost_i,ghost_j,ghost_k)	= mix_mol_mass%cells(i,j,k)
										v_s%cells(ghost_i,ghost_j,ghost_k)				= v_s%cells(i,j,k)
							
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = -v%pr(dim1)%cells(i,j,k)
											else
												if(bc%boundary_types(bound_number)%is_slip()) then
													v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = v%pr(dim1)%cells(i,j,k)
												else
													v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = -v%pr(dim1)%cells(i,j,k)
												end if
											end if
										end do

										do specie_number = 1, this%chem%chem_ptr%species_number
											Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)	= Y%pr(specie_number)%cells(i,j,k)
										end do
			
										if(bc%boundary_types(bound_number)%is_conductive()) then
											wall_temperature = bc%boundary_types(bound_number)%get_wall_temperature()
											T%cells(ghost_i,ghost_j,ghost_k) = wall_temperature
										end if

									case('symmetry_plane')
										! A symmetry plane is a slip reflective boundary: normal velocity is odd,
										! thermodynamic variables, tangential velocities and species are even.
										p%cells(ghost_i,ghost_j,ghost_k)					= p%cells(i,j,k)
										rho%cells(ghost_i,ghost_j,ghost_k)				= rho%cells(i,j,k)
										T%cells(ghost_i,ghost_j,ghost_k)					= T%cells(i,j,k)
										mix_mol_mass%cells(ghost_i,ghost_j,ghost_k)	= mix_mol_mass%cells(i,j,k)
										v_s%cells(ghost_i,ghost_j,ghost_k)				= v_s%cells(i,j,k)
										do dim1 = 1, dimensions
											if(dim1 == dim) then
												v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = -v%pr(dim1)%cells(i,j,k)
											else
												v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = v%pr(dim1)%cells(i,j,k)
											end if
										end do
										do specie_number = 1, this%chem%chem_ptr%species_number
											Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)	= Y%pr(specie_number)%cells(i,j,k)
										end do

									case ('inlet')
										! Prescribed inflow/far-field ghost state for conservative variables.
										farfield_pressure	 = bc%boundary_types(bound_number)%get_farfield_pressure()
										farfield_density	 = bc%boundary_types(bound_number)%get_farfield_density()
										farfield_temperature = bc%boundary_types(bound_number)%get_farfield_temperature()
										farfield_velocity	 = bc%boundary_types(bound_number)%get_farfield_velocity()

										if (farfield_pressure > 0.0_dp) p%cells(ghost_i,ghost_j,ghost_k) = farfield_pressure
										if (farfield_density > 0.0_dp) rho%cells(ghost_i,ghost_j,ghost_k) = farfield_density
										if (farfield_temperature > 0.0_dp) T%cells(ghost_i,ghost_j,ghost_k) = farfield_temperature

										if ((farfield_pressure > 0.0_dp).and. &
											(farfield_density > 0.0_dp).and. &
											(farfield_temperature > 0.0_dp)) then
											mol_mix_from_farfield = farfield_density*farfield_temperature*r_gase_J/farfield_pressure
											mix_mol_mass%cells(ghost_i,ghost_j,ghost_k) = mol_mix_from_farfield
										else
											mix_mol_mass%cells(ghost_i,ghost_j,ghost_k) = mix_mol_mass%cells(i,j,k)
										end if

										if (v_s%cells(ghost_i,ghost_j,ghost_k) <= 0.0_dp) then
											v_s%cells(ghost_i,ghost_j,ghost_k) = v_s%cells(i,j,k)
										end if

										v%pr(dim)%cells(ghost_i,ghost_j,ghost_k) = farfield_velocity
										do dim1 = 1, dimensions
											if (dim1 /= dim) then
												v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = v%pr(dim1)%cells(i,j,k)
											end if
										end do

										spec_summ = 0.0_dp
										do specie_number = 1, this%chem%chem_ptr%species_number
											Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = &
												max(Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k), 0.0_dp)
											spec_summ = spec_summ + Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)
										end do
										if (spec_summ > 0.0_dp) then
											do specie_number = 1, this%chem%chem_ptr%species_number
												Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = &
													Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)/spec_summ
											end do
										else
											do specie_number = 1, this%chem%chem_ptr%species_number
												Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = Y%pr(specie_number)%cells(i,j,k)
											end do
										end if

									case ('outlet')
										! Conservative ghost state for transport/source solvers.
										! True outflow uses zero normal gradients; backflow is prescribed.
										normal_velocity = real(sign, dp)*v%pr(dim)%cells(i,j,k)
										if (normal_velocity >= 0.0_dp) then
											p%cells(ghost_i,ghost_j,ghost_k)					= p%cells(i,j,k)
											rho%cells(ghost_i,ghost_j,ghost_k)				= rho%cells(i,j,k)
											T%cells(ghost_i,ghost_j,ghost_k)					= T%cells(i,j,k)
											mix_mol_mass%cells(ghost_i,ghost_j,ghost_k)	= mix_mol_mass%cells(i,j,k)
											v_s%cells(ghost_i,ghost_j,ghost_k)				= v_s%cells(i,j,k)
											do dim1 = 1, dimensions
												v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = v%pr(dim1)%cells(i,j,k)
											end do
											do specie_number = 1, this%chem%chem_ptr%species_number
												Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = &
													Y%pr(specie_number)%cells(i,j,k)
											end do
										else
											farfield_pressure	 = bc%boundary_types(bound_number)%get_farfield_pressure()
											farfield_density	 = bc%boundary_types(bound_number)%get_farfield_density()
											farfield_temperature = bc%boundary_types(bound_number)%get_farfield_temperature()
											farfield_velocity	 = bc%boundary_types(bound_number)%get_farfield_velocity()
											if (farfield_pressure > 0.0_dp) then
												p%cells(ghost_i,ghost_j,ghost_k) = farfield_pressure
											else
												p%cells(ghost_i,ghost_j,ghost_k) = p%cells(i,j,k)
											end if
											if (farfield_density > 0.0_dp) then
												rho%cells(ghost_i,ghost_j,ghost_k) = farfield_density
											else
												rho%cells(ghost_i,ghost_j,ghost_k) = rho%cells(i,j,k)
											end if
											if (farfield_temperature > 0.0_dp) then
												T%cells(ghost_i,ghost_j,ghost_k) = farfield_temperature
											else
												T%cells(ghost_i,ghost_j,ghost_k) = T%cells(i,j,k)
											end if
											if ((farfield_pressure > 0.0_dp).and. &
												(farfield_density > 0.0_dp).and. &
												(farfield_temperature > 0.0_dp)) then
												mol_mix_from_farfield = farfield_density*farfield_temperature*r_gase_J/farfield_pressure
												mix_mol_mass%cells(ghost_i,ghost_j,ghost_k) = mol_mix_from_farfield
											else
												mix_mol_mass%cells(ghost_i,ghost_j,ghost_k) = mix_mol_mass%cells(i,j,k)
											end if
											if (v_s%cells(ghost_i,ghost_j,ghost_k) <= 0.0_dp) then
												v_s%cells(ghost_i,ghost_j,ghost_k) = v_s%cells(i,j,k)
											end if
											v%pr(dim)%cells(ghost_i,ghost_j,ghost_k) = farfield_velocity
											do dim1 = 1, dimensions
												if (dim1 /= dim) then
													v%pr(dim1)%cells(ghost_i,ghost_j,ghost_k) = v%pr(dim1)%cells(i,j,k)
												end if
											end do
											! Preserve any initialized backflow/far-field composition if present.
											spec_summ = 0.0_dp
											do specie_number = 1, this%chem%chem_ptr%species_number
												Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = &
													max(Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k), 0.0_dp)
												spec_summ = spec_summ + Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)
											end do
											if (spec_summ > 0.0_dp) then
												do specie_number = 1, this%chem%chem_ptr%species_number
													Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = &
														Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k)/spec_summ
												end do
											else
												do specie_number = 1, this%chem%chem_ptr%species_number
													Y%pr(specie_number)%cells(ghost_i,ghost_j,ghost_k) = Y%pr(specie_number)%cells(i,j,k)
												end do
											end if
										end if

									case default
										write(*,*) 'CABARET apply_boundary_conditions_main: unsupported boundary type ', &
											trim(boundary_type_name)
										stop
								end select

							end if
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
	
	!> Compute the global explicit time step from the multidimensional acoustic CFL constraint.
	
	subroutine calculate_time_step(this)

#ifdef mpi
	use MPI
#endif

		class(cabaret_solver)	,intent(inout)	:: this
		
		real(dp)	:: local_min_dt, global_min_dt, delta_t_interm
		real(dp)	:: directional_wave_sum
		integer		:: dimensions
		integer		:: mpi_communicator
		integer		,dimension(3,2)	:: cons_inner_loop
		real(dp)	,dimension(3)	:: cell_size
		integer		:: i,j,k,dim,error

		mpi_communicator	= this%domain%get_mpi_communicator()
		local_min_dt	= huge(1.0_dp)
		global_min_dt	= huge(1.0_dp)

		associate(  v				=> this%v%v_ptr		, &
					v_s				=> this%v_s%s_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)
		
		dimensions			= this%domain%get_domain_dimensions()
		cons_inner_loop		= this%domain%get_local_inner_cells_bounds()
		cell_size			= mesh%get_cell_edges_length()					
					
		! Use a multidimensional explicit CFL estimate.  The old implementation
		! used min(dx)/(|u|+c), which is acceptable on the current uniform grid,
		! but this form is easier to generalize to dx /= dy /= dz and AMR.
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if(bc%bc_markers(i,j,k) == 0) then
				directional_wave_sum = 0.0_dp
				do dim = 1,dimensions
					if (cell_size(dim) > 0.0_dp) then
						directional_wave_sum = directional_wave_sum + (abs(v%pr(dim)%cells(i,j,k)) + v_s%cells(i,j,k)) / cell_size(dim)
					end if
				end do

				if (directional_wave_sum > 0.0_dp) then
					delta_t_interm = 1.0_dp / directional_wave_sum
					local_min_dt = min(local_min_dt, delta_t_interm)
				end if
			end if
		end do
		end do
		end do
	
#ifdef mpi					
		call mpi_allreduce(local_min_dt, global_min_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpi_communicator, error)
#else
		global_min_dt = local_min_dt
#endif

		! If all local cells are inactive, retain the previously configured step.
		if (global_min_dt < huge(1.0_dp)) &
			this%time_step = this%courant_fraction * global_min_dt

		end associate
			
	end subroutine calculate_time_step

	!> Set the user-controlled Courant multiplier.

	subroutine set_CFL_coefficient(this,coefficient)
		class(cabaret_solver)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: coefficient
	
		this%courant_fraction = coefficient
		
	end subroutine
	
	!> Return the current time-step size.
	
	pure function get_time_step(this)
		real(dp)						:: get_time_step
		class(cabaret_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	!> Return the current simulation time.

	pure function get_time(this)
		real(dp)						:: get_time
		class(cabaret_solver)	,intent(in)		:: this

		get_time = this%time
	end function

end module
