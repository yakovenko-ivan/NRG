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

	! Toggle to .true. when debugging conservation/thermodynamic consistency.
	logical, parameter :: cabaret_debug_checks = .false.

	! Optional local Riemann usage at material contacts.  The main production
	! mode in this version is pressure-only projection before entropy locking:
	! use the Riemann solver to repair only p_f when entropy inversion detects
	! an anomalous temperature, leaving CABARET Y_f, s_f and velocities intact.
	logical, parameter :: cabaret_use_riemann_pressure_projection = .false.
	logical, parameter :: cabaret_riemann_pressure_projection_blend = .true.
	logical, parameter :: cabaret_riemann_pressure_projection_print_statistics = .false.
	real(dp), parameter :: cabaret_riemann_pressure_projection_T_ratio = 1.20_dp
	real(dp), parameter :: cabaret_riemann_pressure_projection_T_abs = 100.0_dp
	real(dp), parameter :: cabaret_riemann_pressure_projection_min_improvement = 0.05_dp
	integer, parameter :: cabaret_riemann_pressure_projection_bisection_iterations = 24

	! Keep the older full-state Riemann repair only as an emergency fallback for
	! invalid face states.  Modes: 0 invalid only; 1 all strong material contacts;
	! 2 invalid plus fully supersonic material contacts.  The default here is 0 to
	! avoid making the Riemann solver dominant and smearing the contact.
	logical, parameter :: cabaret_use_riemann_face_repair = .true.
	integer, parameter :: cabaret_riemann_face_repair_mode = 0
	logical, parameter :: cabaret_use_thermally_perfect_riemann = .true.
	logical, parameter :: cabaret_riemann_face_repair_print_statistics = .false.
	real(dp), parameter :: cabaret_riemann_molar_mass_ratio = 1.10_dp
	real(dp), parameter :: cabaret_riemann_supersonic_eps = 1.0e-10_dp

	! Initial CABARET flux variables must be consistent with the cell
	! conservative state.  At discontinuities, especially multicomponent
	! material contacts, use a local 1-D Riemann problem instead of
	! arithmetic averaging.  Smooth faces keep the arithmetic average unless
	! cabaret_riemann_initial_all_interior_faces is enabled for diagnostics.
	logical, parameter :: cabaret_use_riemann_initial_faces = .true.
	logical, parameter :: cabaret_riemann_initial_all_interior_faces = .false.
	real(dp), parameter :: cabaret_initial_riemann_molar_mass_ratio = 1.01_dp
	real(dp), parameter :: cabaret_initial_riemann_pressure_ratio = 1.05_dp
	real(dp), parameter :: cabaret_initial_riemann_velocity_ratio = 0.05_dp
	real(dp), parameter :: cabaret_initial_contact_velocity_eps = 10.0_dp*tiny(1.0_dp)

	! Optional research mode: use an entropy/isentrope proxy K=p/rho**gamma
	! instead of the default linear contact quasi-invariant p-c**2*rho.
	! The default remains the original linear balance-characteristic form.
	logical, parameter :: cabaret_use_entropy_contact_invariant = .false.

	! Source shifts in the characteristic maximum-principle intervals.  Geometry is
	! a true in-step balance source of the local finite-volume gas-dynamic update.
	! The optional effective-source mode evaluates selected physical source solvers
	! once per CABARET time step and interprets their production fields as full-step
	! averaged conservative source rates.  The same rates are used in the predictor,
	! corrector and characteristic limiter shifts; included sources are not applied
	! again in the operator-split source stage.
	logical, parameter :: cabaret_use_geometry_source_shifts = .true.
	logical, parameter :: cabaret_use_effective_physical_sources = .true.

	! Keep conservative effective-source coupling separated from characteristic
	! source shifts. This makes it possible to test the predictor/corrector
	! sequence with all physical g_l shifts disabled.
	logical, parameter :: cabaret_use_effective_physical_source_shifts = .true.
	logical, parameter :: cabaret_include_chemistry_in_gas_step = .true.
	logical, parameter :: cabaret_include_diffusion_in_gas_step = .true.
	logical, parameter :: cabaret_include_heat_transfer_in_gas_step = .true.
	logical, parameter :: cabaret_include_viscosity_in_gas_step = .true.
	logical, parameter :: cabaret_include_particles_in_gas_step = .false.

	! Experimental multicomponent-contact fix. Species and entropy are
	! reconstructed as one material/contact package; the face density is then
	! recovered from p_f, Y_f and the reconstructed mixture entropy instead of
	! from the linear density invariant alone.
	logical, parameter :: cabaret_use_entropy_locked_contact_density = .true.
	logical, parameter :: cabaret_entropy_lock_all_interior_faces = .false.
	real(dp), parameter :: cabaret_entropy_lock_molar_mass_ratio = 1.000001_dp
	real(dp), parameter :: cabaret_entropy_p_ref = 101325.0_dp

	! Very local safety limiter for entropy-locked density closure.  The
	! temperature guard was the first successful diagnostic/fallback, but the
	! production candidate below limits the thermodynamic variable that actually
	! causes the one-cell spike at a material tail:
	!
	!     psi = M(Y)/rho = R_u*T/p .
	!
	! At a nearly isobaric contact, a temperature spike is a spike in psi.  The
	! molar-volume guard clips only the density closure implied by entropy lock;
	! p_f, Y_f, s_f and velocities are left unchanged.
	logical, parameter :: cabaret_use_entropy_lock_temperature_guard = .false.
	logical, parameter :: cabaret_entropy_lock_temperature_guard_print_statistics = .false.
	real(dp), parameter :: cabaret_entropy_lock_temperature_guard_ratio = 1.10_dp
	real(dp), parameter :: cabaret_entropy_lock_temperature_guard_abs = 25.0_dp

	logical, parameter :: cabaret_use_entropy_lock_molar_volume_guard = .true.
	logical, parameter :: cabaret_entropy_lock_molar_volume_guard_print_statistics = .true.

	! Selective activation for the molar-volume guard.  Entropy lock remains active
	! for very weak material jumps, but the guard is applied only to stronger
	! material contacts and only when the entropy-closed state would create a
	! visible high-temperature consequence.  The trigger tolerance is deliberately
	! looser than the final clip tolerance: harmless second-order extrapolations are
	! tolerated, while true one-cell outliers are projected back close to the local
	! two-cell molar-volume interval.
	logical, parameter :: cabaret_entropy_lock_molar_volume_guard_upper_only = .true.
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_molar_mass_ratio = 1.05_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_trigger_rel_tol = 0.05_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_trigger_abs_tol = 0.0_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_clip_rel_tol = 0.01_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_clip_abs_tol = 0.0_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_T_ratio = -1.0_dp
	real(dp), parameter :: cabaret_entropy_lock_molar_volume_guard_T_abs = 0.0_dp

	! Strict thermodynamic inversion bracket used only by entropy->temperature
	! inversion.  It is no longer used as a silent temperature repair/clipping
	! floor or ceiling elsewhere.
	real(dp), parameter :: cabaret_entropy_temperature_bracket_low = 1.0_dp
	real(dp), parameter :: cabaret_entropy_temperature_bracket_high = tables_temperature_ceiling

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

		! Effective conservative source rates used by the optional source-aware
		! CABARET predictor/corrector.  They are evaluated once per full time step
		! before the predictor.  If a source solver internally integrates over dt,
		! the corresponding production field is interpreted as an averaged rate:
		! rho_src=d(rho)/dt, mom_src=d(rho*u)/dt, rhoE_src=d(rho*E)/dt,
		! rhoY_src=d(rho*Y_k)/dt.
		real(dp) ,dimension(:,:,:)	,allocatable    :: rho_src_old, rhoE_src_old
		real(dp)	,dimension(:,:,:,:)	,allocatable	:: mom_src_old, rhoY_src_old

		! Old-time flux-divergence residuals A^n saved during the predictor.
		! In effective-source mode the final conservative update is formed as
		! U^{n+1}=U^n+0.5*dt*(A^n+A^{n+1})+dt*S_eff, rather than by
		! anchoring the corrector on a cached, source-modified half-step state.
		real(dp) ,dimension(:,:,:)	,allocatable    :: rho_rhs_old, rhoE_rhs_old
		real(dp)	,dimension(:,:,:,:)	,allocatable	:: mom_rhs_old, rhoY_rhs_old

		! Flow variables
		real(dp) ,dimension(:,:,:,:,:)	,allocatable    :: v_f, Y_f
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: rho_f, p_f, e_i_f, E_f_f, v_s_f
		! Quasi invariants
		real(dp) ,dimension(:,:,:,:,:)		,allocatable    :: r_inv_corr, q_inv_corr
        real(dp) ,dimension(:,:,:,:,:,:)		,allocatable    :: v_inv_corr

		! Reconstructed material/contact entropy at faces. It is a temporary
		! solver-owned closure variable used to lock rho_f to the same material
		! package as Y_f.
		real(dp) ,dimension(:,:,:,:)		,allocatable    :: s_material_f_new
        
        
	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure	,private	:: apply_boundary_conditions_flow
		procedure				:: solve_problem
		procedure				:: solve_test_problem
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
		procedure	,private	:: initialize_material_entropy_faces
		procedure	,private	:: apply_riemann_pressure_projection
		procedure	,private	:: enforce_material_contact_density_from_entropy
		procedure	,private	:: apply_optional_riemann_face_repair
		procedure	,private	:: update_flow_thermodynamics
		procedure	,private	:: correct_conservative_full_step
		procedure	,private	:: finalize_gas_dynamics_step
		procedure	,private	:: solve_split_physics
		procedure	,private	:: apply_split_sources
		procedure	,private	:: cache_flow_state_for_next_step
		procedure	,private	:: cell_length
		procedure	,private	:: state_is_finite
		procedure	,private	:: check_conservative_state
		procedure	,private	:: check_face_state
	end type

	interface	cabaret_solver_c
		module procedure	constructor
	end interface

contains

	pure real(dp) function effective_gamma_from_state(p_state, rho_state, c_state) result(gamma_eff)
		real(dp), intent(in) :: p_state, rho_state, c_state

		! Strict mode: do not repair pressure/density/gamma. Invalid inputs are
		! allowed to produce IEEE exceptions or propagate NaNs instead of being
		! silently projected to an admissible state.
		gamma_eff = c_state*c_state*rho_state/p_state
	end function effective_gamma_from_state


	pure real(dp) function contact_quasi_invariant(p_state, rho_state, c_state) result(inv)
		real(dp), intent(in) :: p_state, rho_state, c_state
		real(dp) :: gamma_eff

		if (cabaret_use_entropy_contact_invariant) then
			gamma_eff = effective_gamma_from_state(p_state, rho_state, c_state)
			! Strict mode: no p/rho/gamma floors.
			inv = p_state/rho_state**gamma_eff
		else
			! Default linearized balance-characteristic contact quasi-invariant.
			inv = p_state - c_state*c_state*rho_state
		end if
	end function contact_quasi_invariant


	pure real(dp) function density_from_contact_quasi_invariant(p_state, inv, rho_ref, c_ref) result(rho_state)
		real(dp), intent(in) :: p_state, inv, rho_ref, c_ref
		real(dp) :: gamma_eff

		if (cabaret_use_entropy_contact_invariant) then
			gamma_eff = effective_gamma_from_state(p_state, rho_ref, c_ref)
			rho_state = (p_state/inv)**(1.0_dp/gamma_eff)
		else
			rho_state = (p_state - inv)/(c_ref*c_ref)
		end if
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
		! predictor/corrector update.  In strict/FPE mode, zero denominators are
		! not hidden by a floor: they should stop the run at the offending state.
		coeff = 2.0_dp*real(nu - 1, dp)/(r_plus**(nu - 1) + r_minus**(nu - 1)) * &
			(r_plus**(nu - 1) - r_minus**(nu - 1))/(r_plus - r_minus)
	end function finite_volume_geometry_coefficient


	pure real(dp) function limit_quasi_invariant(value, left_value, half_value, right_value, source_shift, alpha, lower_clip, upper_clip) result(limited)
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
		real(dp)				:: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, v_l, v_r
		real(dp)				:: c_l, c_r, p_ratio, velocity_scale
		real(dp)				:: p_floor, rho_floor, u_face_init, u_contact_init
		real(dp)	,dimension(:)	,allocatable :: Y_left_riemann, Y_right_riemann, Y_face_riemann
		logical				:: use_left_contact_state, use_riemann_initial_face
		integer	,dimension(3,2)	:: cons_allocation_bounds, flow_allocation_bounds
		integer	,dimension(3,2)	:: flow_inner_loop, loop
		
		integer					:: dimensions, species_number
		integer					:: i, j, k, dim, dim1, spec
		real(dp)	,dimension(3)	:: cell_size

		cons_allocation_bounds		= manager%domain%get_local_utter_cells_bounds()
		flow_allocation_bounds		= manager%domain%get_local_utter_faces_bounds()
		dimensions					= manager%domain%get_domain_dimensions()

		species_number			= manager%chemistry%chem_ptr%species_number
		allocate(Y_left_riemann(species_number), Y_right_riemann(species_number), Y_face_riemann(species_number))
		
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
        
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'density')
		constructor%rho%s_ptr				=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'temperature')
		constructor%T%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'pressure')
		constructor%p%s_ptr					=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'full_energy')
		constructor%E_f%s_ptr				=> scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'internal_energy')
		constructor%e_i%s_ptr				=> scal_c_ptr%s_ptr		
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_c_ptr%s_ptr		
		
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
			constructor%diff_solver			= diffusion_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_diffusion')
			constructor%E_f_prod_diff%s_ptr			=> scal_c_ptr%s_ptr			
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_c_ptr%v_ptr
		end if

		if (constructor%heat_trans_flag) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr			=> scal_c_ptr%s_ptr
		end if

		if(constructor%viscosity_flag) then
			constructor%viscosity_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_viscosity')
			constructor%E_f_prod_visc%s_ptr			=> scal_c_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_c_ptr%v_ptr
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
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%rho_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr                
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter						
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)						!# Continuum particles solver								
				constructor%v_prod_particles(particles_phase_counter)%v_ptr		=> vect_c_ptr%v_ptr						!# Continuum particles solver
				write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%Y_prod_particles(particles_phase_counter)%v_ptr		=> vect_c_ptr%v_ptr                
			end do		
		end if
	
		call manager%get_flow_field_pointer_by_name(scal_f_ptr,vect_f_ptr,tens_f_ptr,'adiabatic_index_flow')
		constructor%gamma_f_new%s_ptr	=> scal_f_ptr%s_ptr				
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

		allocate(constructor%s_material_f_new(	dimensions						, &
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
			constructor%s_material_f_new(:,:,:,:) = 0.0_dp
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
								gamma_l = constructor%gamma%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								gamma_r = constructor%gamma%s_ptr%cells(i,j,k)
								v_l     = constructor%v%v_ptr%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								v_r     = constructor%v%v_ptr%pr(dim)%cells(i,j,k)
								c_l     = constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								c_r     = constructor%v_s%s_ptr%cells(i,j,k)

								p_ratio = max(p_l,p_r) / min(p_l,p_r)
								velocity_scale = 0.5_dp*(abs(c_l) + abs(c_r))
								use_riemann_initial_face = cabaret_use_riemann_initial_faces .and. &
									(cabaret_riemann_initial_all_interior_faces .or. &
									 molar_mass_ratio > cabaret_initial_riemann_molar_mass_ratio .or. &
									 p_ratio > cabaret_initial_riemann_pressure_ratio .or. &
									 abs(v_l - v_r) > cabaret_initial_riemann_velocity_ratio*velocity_scale)

								if (use_riemann_initial_face) then
									p_floor = 1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r))
									rho_floor = 1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r))
									do spec = 1, species_number
										Y_left_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
										Y_right_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
									end do

									if (cabaret_use_thermally_perfect_riemann) then
										call initial_riemann%set_thermally_perfect_parameters(constructor%thermo%thermo_ptr, &
											rho_l, rho_r, p_l, p_r, v_l, v_r, Y_left_riemann, Y_right_riemann, p_floor, rho_floor)
									else
										call initial_riemann%set_parameters(rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, v_l, v_r, p_floor, rho_floor)
									end if
									call initial_riemann%solve()

									if (initial_riemann%get_success()) then
										u_face_init = initial_riemann%get_velocity()
										u_contact_init = initial_riemann%get_contact_velocity()
										use_left_contact_state = (u_contact_init >= 0.0_dp)
										if (cabaret_use_thermally_perfect_riemann) call initial_riemann%get_mass_fractions(Y_face_riemann)

										constructor%p_f(dim,i,j,k)   = initial_riemann%get_pressure()
										constructor%rho_f(dim,i,j,k) = initial_riemann%get_density()
										constructor%v_f(dim,dim,i,j,k) = u_face_init
									else
										! Robust fallback: keep a single material side rather than an
										! averaged composition.  At an initially motionless diaphragm,
										! the pressure jump gives the first contact direction.
										u_face_init = 0.5_dp*(v_l + v_r)
										if (abs(u_face_init) > cabaret_initial_contact_velocity_eps) then
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
											if (initial_riemann%get_success() .and. cabaret_use_thermally_perfect_riemann) then
												constructor%Y_f(spec,dim,i,j,k) = Y_face_riemann(spec)
											else
												constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
											end if
										end do
										do dim1 = 1, dimensions
										if (dim1 /= dim) then
											constructor%v_f(dim1,dim,i,j,k) = constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
										end if
										end do
									else
										constructor%E_f_f(dim,i,j,k) = constructor%E_f%s_ptr%cells(i,j,k)
										constructor%v_s_f(dim,i,j,k) = constructor%v_s%s_ptr%cells(i,j,k)
										do spec = 1, species_number
											if (initial_riemann%get_success() .and. cabaret_use_thermally_perfect_riemann) then
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
		constructor%initial_time_step   =   manager%solver_options%get_initial_time_step()
		
		call constructor%state_eq%apply_state_equation_flow_variables_for_IC()

		call constructor%apply_boundary_conditions_main()
		
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

		call cabaret_timer%tic()

		! CABARET gas-dynamics stage: predictor -> face reconstruction -> corrector.
		! In effective-source mode selected physical source solvers are evaluated
		! once here.  Their full-step production fields are used as averaged
		! conservative source rates in both CABARET half-steps and in the
		! characteristic limiter shifts.
		call this%prepare_gas_dynamics_step()
		call this%evaluate_effective_physical_source_rates()

		call cabaret_gas_dynamics_timer%tic()
		call this%predict_conservative_half_step()
		call cabaret_gas_dynamics_timer%toc()

		call this%update_cell_thermodynamics()
		call this%check_conservative_state('after predictor')

		call cabaret_gas_dynamics_timer%tic()
		call this%reconstruct_acoustic_face_state()
		call this%reconstruct_contact_family_face_state()
		call this%finish_face_reconstruction()
		call this%apply_riemann_pressure_projection()
		call this%enforce_material_contact_density_from_entropy()
		call this%apply_optional_riemann_face_repair()
		call cabaret_gas_dynamics_timer%toc()

		call this%update_flow_thermodynamics()
		call this%check_face_state('after face EOS')

		call cabaret_gas_dynamics_timer%tic()
		call this%correct_conservative_full_step()
		call this%finalize_gas_dynamics_step()
		call cabaret_gas_dynamics_timer%toc()

		! Operator-split physics: heat transfer, diffusion, viscosity, chemistry, particles.
		call this%solve_split_physics()

		call cabaret_gas_dynamics_timer%tic()
		call this%apply_split_sources()
		call cabaret_gas_dynamics_timer%toc()

		call this%update_cell_thermodynamics()
		call this%check_conservative_state('after split sources')

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

	subroutine prepare_gas_dynamics_step(this)
		class(cabaret_solver), intent(inout) :: this

		call this%apply_boundary_conditions_main()
		call this%exchange_conservative_state()
		call this%cache_conservative_state()
	end subroutine prepare_gas_dynamics_step

	subroutine evaluate_effective_physical_source_rates(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, spec
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

		if (.not. cabaret_use_effective_physical_sources) return

		! Evaluate selected physical source solvers once per full CABARET time step.
		! Existing production fields are interpreted as effective conservative source
		! rates, or as averaged rates if the solver internally integrates over dt.
		if (this%heat_trans_flag .and. cabaret_include_heat_transfer_in_gas_step) then
			call cabaret_heattransfer_timer%tic()
			call this%heat_trans_solver%solve_heat_transfer(this%time_step)
			call cabaret_heattransfer_timer%toc(new_iter=.true.)
		end if

		if (this%diffusion_flag .and. cabaret_include_diffusion_in_gas_step) then
			call cabaret_diffusion_timer%tic()
			call this%diff_solver%solve_diffusion(this%time_step)
			call cabaret_diffusion_timer%toc(new_iter=.true.)
		end if

		if (this%viscosity_flag .and. cabaret_include_viscosity_in_gas_step) then
			call cabaret_viscosity_timer%tic()
			call this%viscosity_solver%solve_viscosity(this%time_step)
			call cabaret_viscosity_timer%toc(new_iter=.true.)
		end if

		if (this%reactive_flag .and. cabaret_include_chemistry_in_gas_step) then
			call cabaret_chemistry_timer%tic()
			call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
			call cabaret_chemistry_timer%toc(new_iter=.true.)
		end if

		associate(  E_f_prod_chem  => this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat  => this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc  => this%E_f_prod_visc%s_ptr	, &
					E_f_prod_diff  => this%E_f_prod_diff%s_ptr	, &
					Y_prod_chem    => this%Y_prod_chem%v_ptr		, &
					Y_prod_diff    => this%Y_prod_diff%v_ptr		, &
					v_prod_visc    => this%v_prod_visc%v_ptr		, &
					bc             => this%boundary%bc_ptr)

		!$omp parallel default(shared) private(i,j,k,dim,spec)
		!$omp do collapse(3) schedule(guided)
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			if (bc%bc_markers(i,j,k) == 0) then
				if (this%reactive_flag .and. cabaret_include_chemistry_in_gas_step) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_chem%cells(i,j,k)
					do spec = 1, species_number
						rhoY_src(spec,i,j,k) = rhoY_src(spec,i,j,k) + Y_prod_chem%pr(spec)%cells(i,j,k)
					end do
				end if

				if (this%heat_trans_flag .and. cabaret_include_heat_transfer_in_gas_step) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_heat%cells(i,j,k)
				end if

				if (this%diffusion_flag .and. cabaret_include_diffusion_in_gas_step) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_diff%cells(i,j,k)
					do spec = 1, species_number
						rhoY_src(spec,i,j,k) = rhoY_src(spec,i,j,k) + Y_prod_diff%pr(spec)%cells(i,j,k)
					end do
				end if

				if (this%viscosity_flag .and. cabaret_include_viscosity_in_gas_step) then
					rhoE_src(i,j,k) = rhoE_src(i,j,k) + E_f_prod_visc%cells(i,j,k)
					do dim = 1, dimensions
						mom_src(dim,i,j,k) = mom_src(dim,i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k)
					end do
				end if
			end if
		end do
		end do
		end do
		!$omp end do nowait
		!$omp end parallel
		end associate

		end associate

	end subroutine evaluate_effective_physical_source_rates


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


	subroutine physical_source_contact_shifts(this, rho_state, p_state, c_state, E_state, vel_state, Y_state, &
			S_rho, S_mom, S_rhoE, S_rhoY, g_Y, g_entropy)
		class(cabaret_solver), intent(in) :: this
		real(dp), intent(in) :: rho_state, p_state, c_state, E_state
		real(dp), dimension(:), intent(in) :: vel_state, Y_state, S_mom, S_rhoY
		real(dp), intent(in) :: S_rho, S_rhoE
		real(dp), dimension(:), intent(out) :: g_Y
		real(dp), intent(out) :: g_entropy

		integer :: spec, species_number
		real(dp) :: f_e, f_p, source_norm, probe_dt
		real(dp) :: rho_trial, p_trial, T_loc, T_trial, s_loc, s_trial
		real(dp), dimension(size(Y_state)) :: Y_trial

		species_number = size(Y_state)

		call this%thermo%thermo_ptr%pressure_source_rate_from_conservative_sources( &
			rho_state, p_state, E_state, vel_state, Y_state, S_rho, S_mom, S_rhoE, S_rhoY, &
			f_p, f_e, g_Y)

		! f_p already contains the multicomponent dp/dY contribution for sources that change composition.
		source_norm = abs(S_rho)/rho_state + abs(f_p)/p_state
		do spec = 1, species_number
			source_norm = source_norm + abs(g_Y(spec))
		end do

		if (source_norm == 0.0_dp) then
			g_entropy = 0.0_dp
			return
		end if

		probe_dt = min(0.5_dp*this%time_step, 1.0e-8_dp/source_norm)
		rho_trial = rho_state + probe_dt*S_rho
		p_trial   = p_state   + probe_dt*f_p
		do spec = 1, species_number
			Y_trial(spec) = Y_state(spec) + probe_dt*g_Y(spec)
		end do

		! Do not clip the source probe.  If the effective source drives the local
		! thermodynamic perturbation outside the admissible domain, stop at the cause
		! rather than silently estimating entropy from a repaired state.
		if (rho_trial <= 0.0_dp) error stop 'Effective source entropy probe produced non-positive density'
		if (p_trial <= 0.0_dp) error stop 'Effective source entropy probe produced non-positive pressure'
		do spec = 1, species_number
			if (Y_trial(spec) < 0.0_dp) error stop 'Effective source entropy probe produced negative mass fraction'
		end do

		T_loc   = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_state, rho_state, Y_state)
		T_trial = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_trial, rho_trial, Y_trial)
		s_loc   = this%thermo%thermo_ptr%mixture_specific_entropy(T_loc, p_state, Y_state, cabaret_entropy_p_ref)
		s_trial = this%thermo%thermo_ptr%mixture_specific_entropy(T_trial, p_trial, Y_trial, cabaret_entropy_p_ref)
		g_entropy = (s_trial - s_loc)/probe_dt

	end subroutine physical_source_contact_shifts


	subroutine update_cell_thermodynamics(this)
		class(cabaret_solver), intent(inout) :: this

		call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation()
		call cabaret_eos_timer%toc(new_iter=.true.)
	end subroutine update_cell_thermodynamics

	subroutine update_flow_thermodynamics(this)
		class(cabaret_solver), intent(inout) :: this

		call cabaret_eos_timer%tic()
		call this%state_eq%apply_state_equation_flow_variables()
		call cabaret_eos_timer%toc(new_iter=.true.)
	end subroutine update_flow_thermodynamics

	subroutine finish_face_reconstruction(this)
		class(cabaret_solver), intent(inout) :: this

		! For the original CABARET half-step corrector, rho_old/v_old/E_f_old/Y_old
		! are deliberately refreshed to the predicted n+1/2 state.  In the
		! effective-source full-corrector variant they must remain equal to the
		! beginning-of-step conservative state U^n, because the final update is
		! assembled as U^n+0.5*dt*(A^n+A^{n+1})+dt*S_eff.
		if (.not. cabaret_use_effective_physical_sources) then
			call this%cache_conservative_state()
		end if
		call this%normalize_face_mass_fractions()
		call this%exchange_face_primitive_state()
	end subroutine finish_face_reconstruction



	subroutine apply_riemann_pressure_projection(this)
		class(cabaret_solver), intent(inout) :: this

		type(riemann_solver) :: riemann
		integer :: dimensions, species_number
		integer :: i, j, k, dim, dim1, spec, iter
		integer :: il, jl, kl
		integer :: attempt_count, success_count, failure_count, rejected_count
		integer, dimension(3,2) :: flow_inner_loop, loop
		real(dp) :: mol_l, mol_r, mol_ratio
		real(dp) :: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r
		real(dp) :: v_l, v_r, p_floor, rho_floor
		real(dp) :: p_cab, rho_cab, p_riemann, p_trial, p_new
		real(dp) :: log_p_cab, log_p_riemann, theta_low, theta_high, theta_mid
		real(dp) :: T_guess, T_state, T_limit, T_riemann, T_trial, T_new
		real(dp) :: old_error, new_error, improvement_limit
		real(dp) :: y_sum
		real(dp), dimension(:), allocatable :: Y_left, Y_right, Y_face
		logical :: invalid_face, contact_trigger, temperature_trigger, accept_pressure

		if (.not. cabaret_use_riemann_pressure_projection) return
		if (.not. cabaret_use_entropy_locked_contact_density) return
		if (.not. allocated(this%s_material_f_new)) return

		dimensions = this%domain%get_domain_dimensions()
		species_number = this%chem%chem_ptr%species_number
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()
		attempt_count = 0
		success_count = 0
		failure_count = 0
		rejected_count = 0

		associate( rho       => this%rho%s_ptr       , &
					p         => this%p%s_ptr         , &
					gamma     => this%gamma%s_ptr     , &
					v         => this%v%v_ptr         , &
					Y         => this%Y%v_ptr         , &
					rho_f_new => this%rho_f_new%s_ptr , &
					p_f_new   => this%p_f_new%s_ptr   , &
					Y_f_new   => this%Y_f_new%v_ptr   , &
					s_f       => this%s_material_f_new , &
					bc        => this%boundary%bc_ptr )

	!$omp parallel default(shared) private(riemann,i,j,k,dim,dim1,spec,iter,il,jl,kl,loop, &
	!$omp& mol_l,mol_r,mol_ratio,rho_l,rho_r,p_l,p_r,gamma_l,gamma_r,v_l,v_r,p_floor,rho_floor, &
	!$omp& p_cab,rho_cab,p_riemann,p_trial,p_new,log_p_cab,log_p_riemann,theta_low,theta_high, &
	!$omp& theta_mid,T_guess,T_state,T_limit,T_riemann,T_trial,T_new,old_error,new_error, &
	!$omp& improvement_limit,y_sum,Y_left,Y_right,Y_face,invalid_face,contact_trigger, &
	!$omp& temperature_trigger,accept_pressure) &
	!$omp& reduction(+:attempt_count,success_count,failure_count,rejected_count)
		call riemann%clear(reset_counter=.true.)
		allocate(Y_left(species_number), Y_right(species_number), Y_face(species_number))

		do dim = 1, dimensions
			loop = flow_inner_loop
			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

		!$omp do collapse(3) schedule(static)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				il = i - I_m(dim,1)
				jl = j - I_m(dim,2)
				kl = k - I_m(dim,3)

				if (bc%bc_markers(i,j,k) /= 0 .or. bc%bc_markers(il,jl,kl) /= 0) cycle

				mol_l = 0.0_dp
				mol_r = 0.0_dp
				y_sum = 0.0_dp
				do spec = 1, species_number
					Y_left(spec) = Y%pr(spec)%cells(il,jl,kl)
					Y_right(spec) = Y%pr(spec)%cells(i,j,k)
					Y_face(spec) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dp)
					y_sum = y_sum + Y_face(spec)
					mol_l = mol_l + Y_left(spec) / this%thermo%thermo_ptr%molar_masses(spec)
					mol_r = mol_r + Y_right(spec) / this%thermo%thermo_ptr%molar_masses(spec)
				end do
				if (y_sum > 0.0_dp) Y_face = Y_face/y_sum

				mol_l = 1.0_dp / mol_l
				mol_r = 1.0_dp / mol_r
				mol_ratio = max(mol_l, mol_r) / min(mol_l, mol_r)
				contact_trigger = (mol_ratio > cabaret_riemann_molar_mass_ratio)
				if (.not. contact_trigger) cycle

				p_cab = p_f_new%cells(dim,i,j,k)
				rho_cab = rho_f_new%cells(dim,i,j,k)
				invalid_face = (p_cab <= 0.0_dp) .or. (rho_cab <= 0.0_dp)

				if (invalid_face) then
					write(*,*) 'CABARET strict pressure projection: invalid face state at ', dim, i, j, k
					write(*,*) '  p_cab, rho_cab = ', p_cab, rho_cab
					error stop 'CABARET strict pressure projection: invalid face state'
				end if

				T_guess = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_cab, rho_cab, Y_face)
				T_state = this%thermo%thermo_ptr%temperature_from_entropy_pressure_Y(s_f(dim,i,j,k), p_cab, Y_face, T_guess, &
			cabaret_entropy_temperature_bracket_low, cabaret_entropy_temperature_bracket_high, cabaret_entropy_p_ref)

				T_limit = max(T_guess*cabaret_riemann_pressure_projection_T_ratio, &
					T_guess + cabaret_riemann_pressure_projection_T_abs)
				temperature_trigger = (T_state > T_limit)
				if (.not. temperature_trigger) cycle

				rho_l = rho%cells(il,jl,kl)
				rho_r = rho%cells(i,j,k)
				p_l = p%cells(il,jl,kl)
				p_r = p%cells(i,j,k)
				gamma_l = gamma%cells(il,jl,kl)
				gamma_r = gamma%cells(i,j,k)
				v_l = v%pr(dim)%cells(il,jl,kl)
				v_r = v%pr(dim)%cells(i,j,k)

				p_floor = 1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r))
				rho_floor = 1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r))

				attempt_count = attempt_count + 1
				if (cabaret_use_thermally_perfect_riemann) then
					call riemann%set_thermally_perfect_parameters(this%thermo%thermo_ptr, &
						rho_l, rho_r, p_l, p_r, v_l, v_r, Y_left, Y_right, p_floor, rho_floor)
				else
					call riemann%set_parameters(rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, v_l, v_r, p_floor, rho_floor)
				end if
				call riemann%solve()

				if (.not. riemann%get_success()) then
					failure_count = failure_count + 1
					cycle
				end if

				p_riemann = riemann%get_pressure()
				if (p_riemann <= 0.0_dp) then
					write(*,*) 'CABARET strict pressure projection: non-positive Riemann pressure at ', dim, i, j, k
					write(*,*) '  p_riemann = ', p_riemann
					error stop 'CABARET strict pressure projection: non-positive Riemann pressure'
				end if

				T_riemann = this%thermo%thermo_ptr%temperature_from_entropy_pressure_Y(s_f(dim,i,j,k), p_riemann, Y_face, T_guess, &
			cabaret_entropy_temperature_bracket_low, cabaret_entropy_temperature_bracket_high, cabaret_entropy_p_ref)
				if (T_riemann >= T_state) then
					rejected_count = rejected_count + 1
					cycle
				end if

				p_new = p_riemann
				T_new = T_riemann
				if (cabaret_riemann_pressure_projection_blend .and. (.not. invalid_face) .and. (T_riemann <= T_limit)) then
					log_p_cab = log(p_cab)
					log_p_riemann = log(p_riemann)
					theta_low = 0.0_dp
					theta_high = 1.0_dp
					do iter = 1, cabaret_riemann_pressure_projection_bisection_iterations
						theta_mid = 0.5_dp*(theta_low + theta_high)
						p_trial = exp((1.0_dp-theta_mid)*log_p_cab + theta_mid*log_p_riemann)
						T_trial = this%thermo%thermo_ptr%temperature_from_entropy_pressure_Y(s_f(dim,i,j,k), p_trial, Y_face, T_guess, &
			cabaret_entropy_temperature_bracket_low, cabaret_entropy_temperature_bracket_high, cabaret_entropy_p_ref)
						if (T_trial <= T_limit) then
							theta_high = theta_mid
						else
							theta_low = theta_mid
						end if
					end do
					p_new = exp((1.0_dp-theta_high)*log_p_cab + theta_high*log_p_riemann)
					T_new = this%thermo%thermo_ptr%temperature_from_entropy_pressure_Y(s_f(dim,i,j,k), p_new, Y_face, T_guess, &
			cabaret_entropy_temperature_bracket_low, cabaret_entropy_temperature_bracket_high, cabaret_entropy_p_ref)
				end if

				old_error = abs(T_state - T_guess)
				new_error = abs(T_new - T_guess)
				improvement_limit = (1.0_dp - cabaret_riemann_pressure_projection_min_improvement)*old_error
				accept_pressure = (p_new > 0.0_dp) .and. &
					((new_error < improvement_limit) .or. (T_new <= T_limit))

				if (accept_pressure) then
					p_f_new%cells(dim,i,j,k) = p_new
					success_count = success_count + 1
				else
					rejected_count = rejected_count + 1
				end if
			end do
			end do
			end do
		!$omp end do nowait
		end do

		deallocate(Y_left, Y_right, Y_face)
	!$omp end parallel

		end associate

		if (cabaret_riemann_pressure_projection_print_statistics .and. attempt_count > 0) then
			print *, 'CABARET Riemann pressure projection: attempts/success/failure/rejected = ', &
				attempt_count, success_count, failure_count, rejected_count
		end if

		call this%exchange_face_pressure_density()
	end subroutine apply_riemann_pressure_projection

	subroutine apply_optional_riemann_face_repair(this)
		class(cabaret_solver), intent(inout) :: this

		type(riemann_solver) :: riemann
		integer :: dimensions, species_number
		integer :: i, j, k, dim, dim1, spec
		integer :: il, jl, kl
		integer :: repair_count, success_count, failure_count
		integer, dimension(3,2) :: flow_inner_loop, loop
		real(dp) :: mol_l, mol_r, mol_ratio
		real(dp) :: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r
		real(dp) :: v_l, v_r, c_l, c_r, mach_margin
		real(dp) :: p_floor, rho_floor, u_star, u_contact, u_selector
		real(dp) :: T_sample, y_sum
		real(dp), dimension(:), allocatable :: Y_left, Y_right, Y_sample
		logical :: invalid_face, contact_trigger, repair_trigger, use_left
		logical :: supersonic_positive, supersonic_negative

		if (.not. cabaret_use_riemann_face_repair) return

		dimensions = this%domain%get_domain_dimensions()
		species_number = this%chem%chem_ptr%species_number
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()
		repair_count = 0
		success_count = 0
		failure_count = 0

		associate( rho       => this%rho%s_ptr       , &
					p         => this%p%s_ptr         , &
					gamma     => this%gamma%s_ptr     , &
					v_s       => this%v_s%s_ptr       , &
					v         => this%v%v_ptr         , &
					Y         => this%Y%v_ptr         , &
					rho_f_new => this%rho_f_new%s_ptr , &
					p_f_new   => this%p_f_new%s_ptr   , &
					v_f_new   => this%v_f_new%v_ptr   , &
					Y_f_new   => this%Y_f_new%v_ptr   , &
					bc        => this%boundary%bc_ptr )

	!$omp parallel default(shared) private(riemann,i,j,k,dim,dim1,spec,il,jl,kl,loop,mol_l,mol_r,mol_ratio, &
	!$omp& rho_l,rho_r,p_l,p_r,gamma_l,gamma_r,v_l,v_r,c_l,c_r,mach_margin,p_floor,rho_floor, &
	!$omp& u_star,u_contact,u_selector,T_sample,y_sum,Y_left,Y_right,Y_sample,invalid_face,contact_trigger, &
	!$omp& repair_trigger,use_left,supersonic_positive,supersonic_negative) &
	!$omp& reduction(+:repair_count,success_count,failure_count)
		call riemann%clear(reset_counter=.true.)
		allocate(Y_left(species_number), Y_right(species_number), Y_sample(species_number))

		do dim = 1, dimensions
			loop = flow_inner_loop
			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

		!$omp do collapse(3) schedule(static)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				il = i - I_m(dim,1)
				jl = j - I_m(dim,2)
				kl = k - I_m(dim,3)

				! This repair is for interior faces only. Physical boundaries are
				! handled by apply_boundary_conditions_flow().
				if (bc%bc_markers(i,j,k) /= 0 .or. bc%bc_markers(il,jl,kl) /= 0) cycle

				mol_l = 0.0_dp
				mol_r = 0.0_dp
				do spec = 1, species_number
					Y_left(spec) = Y%pr(spec)%cells(il,jl,kl)
					Y_right(spec) = Y%pr(spec)%cells(i,j,k)
					mol_l = mol_l + Y_left(spec) / this%thermo%thermo_ptr%molar_masses(spec)
					mol_r = mol_r + Y_right(spec) / this%thermo%thermo_ptr%molar_masses(spec)
				end do
				mol_l = 1.0_dp / mol_l
				mol_r = 1.0_dp / mol_r
				mol_ratio = max(mol_l, mol_r) / min(mol_l, mol_r)

				rho_l = rho%cells(il,jl,kl)
				rho_r = rho%cells(i,j,k)
				p_l = p%cells(il,jl,kl)
				p_r = p%cells(i,j,k)
				gamma_l = gamma%cells(il,jl,kl)
				gamma_r = gamma%cells(i,j,k)
				v_l = v%pr(dim)%cells(il,jl,kl)
				v_r = v%pr(dim)%cells(i,j,k)
				c_l = v_s%cells(il,jl,kl)
				c_r = v_s%cells(i,j,k)

				invalid_face = (p_f_new%cells(dim,i,j,k) <= 0.0_dp) .or. &
							   (rho_f_new%cells(dim,i,j,k) <= 0.0_dp)
				contact_trigger = (mol_ratio > cabaret_riemann_molar_mass_ratio)

				mach_margin = cabaret_riemann_supersonic_eps*max(c_l, c_r)
				supersonic_positive = (v_l - c_l > mach_margin) .and. (v_r - c_r > mach_margin)
				supersonic_negative = (v_l + c_l < -mach_margin) .and. (v_r + c_r < -mach_margin)

				select case (cabaret_riemann_face_repair_mode)
				case (0)
					repair_trigger = invalid_face
				case (1)
					repair_trigger = invalid_face .or. contact_trigger
				case default
					repair_trigger = invalid_face .or. &
						(contact_trigger .and. (supersonic_positive .or. supersonic_negative))
				end select

				if (.not. repair_trigger) cycle
				repair_count = repair_count + 1

				p_floor = 1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r))
				rho_floor = 1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r))

				if (cabaret_use_thermally_perfect_riemann) then
					call riemann%set_thermally_perfect_parameters(this%thermo%thermo_ptr, &
						rho_l, rho_r, p_l, p_r, v_l, v_r, Y_left, Y_right, p_floor, rho_floor)
				else
					call riemann%set_parameters(rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, v_l, v_r, p_floor, rho_floor)
				end if
				call riemann%solve()

				if (riemann%get_success()) then
					success_count = success_count + 1
					! get_velocity() is the state sampled at x/t=0.  The contact
					! speed selects the tangential velocity/material side when the
					! legacy gamma-law API is used.
					u_star = riemann%get_velocity()
					u_contact = riemann%get_contact_velocity()
					use_left = (u_contact >= 0.0_dp)
					rho_f_new%cells(dim,i,j,k) = riemann%get_density()
					p_f_new%cells(dim,i,j,k) = riemann%get_pressure()
					v_f_new%pr(dim)%cells(dim,i,j,k) = u_star

					if (cabaret_use_thermally_perfect_riemann) then
						call riemann%get_mass_fractions(Y_sample)
						y_sum = sum(max(Y_sample,0.0_dp))
						if (y_sum <= 0.0_dp) then
							if (use_left) then
								Y_sample = Y_left
							else
								Y_sample = Y_right
							end if
						end if
					else
						if (use_left) then
							Y_sample = Y_left
						else
							Y_sample = Y_right
						end if
					end if
				else
					failure_count = failure_count + 1
					if (invalid_face) then
						write(*,*) 'CABARET strict Riemann face repair: invalid face and Riemann solver failed at ', dim, i, j, k
						write(*,*) '  p_f_new, rho_f_new = ', p_f_new%cells(dim,i,j,k), rho_f_new%cells(dim,i,j,k)
						error stop 'CABARET strict Riemann face repair: invalid face and Riemann solver failed'
					end if
					cycle
				end if

				do spec = 1, species_number
					Y_f_new%pr(spec)%cells(dim,i,j,k) = Y_sample(spec)
				end do

				if (use_left) then
					do dim1 = 1, dimensions
						if (dim1 /= dim) v_f_new%pr(dim1)%cells(dim,i,j,k) = v%pr(dim1)%cells(il,jl,kl)
					end do
				else
					do dim1 = 1, dimensions
						if (dim1 /= dim) v_f_new%pr(dim1)%cells(dim,i,j,k) = v%pr(dim1)%cells(i,j,k)
					end do
				end if

				if (allocated(this%s_material_f_new)) then
					if (riemann%get_success() .and. cabaret_use_thermally_perfect_riemann .and. &
						(riemann%get_temperature() > 0.0_dp)) then
						T_sample = riemann%get_temperature()
					else
						T_sample = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_f_new%cells(dim,i,j,k), &
							rho_f_new%cells(dim,i,j,k), Y_sample)
					end if
					this%s_material_f_new(dim,i,j,k) = this%thermo%thermo_ptr%mixture_specific_entropy(T_sample, &
						p_f_new%cells(dim,i,j,k), Y_sample, cabaret_entropy_p_ref)
				end if
			end do
			end do
			end do
		!$omp end do nowait
		end do

		deallocate(Y_left, Y_right, Y_sample)
	!$omp end parallel

        end associate

		if (cabaret_riemann_face_repair_print_statistics .and. repair_count > 0) then
			print *, 'CABARET Riemann face repair: attempts/success/failure = ', &
				repair_count, success_count, failure_count
		end if

		call this%normalize_face_mass_fractions()
		call this%exchange_face_primitive_state()
	end subroutine apply_optional_riemann_face_repair

	subroutine predict_conservative_half_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, spec
		integer, dimension(3,2) :: cons_inner_loop
		real(dp), dimension(3) :: cell_size
		character(len=20) :: coordinate_system
		real(dp) :: spec_summ, mean_higher, mean_lower, r
		real(dp), dimension(:), allocatable :: rhoY_new
		real(dp) :: Max_v_s, Min_v_s

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		coordinate_system = this%domain%get_coordinate_system_name()
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()
		Max_v_s = -huge(1.0_dp)
		Min_v_s =  huge(1.0_dp)

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

        !$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r,rhoY_new)
		allocate(rhoY_new(species_number))

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
						mean_higher	=  rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * y_f(spec,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f(dim,i,j,k) * y_f(spec,dim,i,j,k) * v_f(dim,dim,i,j,k)
						rhoY_new(spec)	=  rhoY_new(spec) - (mean_higher - mean_lower ) /cell_size(1)
						if(dim == 1) then
							rhoY_new(spec)	=  rhoY_new(spec)   - 2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
															    * 0.5_dp * (mean_higher + mean_lower)														&
															    * ((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))											&
															    / ((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1)))
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
					
					mom_rhs_old(dim,i,j,k) = v%pr(dim)%cells(i,j,k)
					v%pr(dim)%cells(i,j,k)	=	rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dp*this%time_step * &
						(mom_rhs_old(dim,i,j,k) + mom_src_old(dim,i,j,k)) 
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


	subroutine reconstruct_acoustic_face_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, dim2, plus, sign, bound_number, spec
		integer, dimension(3,2) :: cons_inner_loop, cons_utter_loop, flow_inner_loop, loop
		real(dp), dimension(3) :: cell_size, characteristic_speed
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
					
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,dim2,plus,loop,G_half,G_half_old,G_half_lower,G_half_higher,r_inv,R_inv_half,R_inv_old,q_inv,Q_inv_half,Q_inv_old,v_inv,v_inv_half,v_inv_old,r_inv_new,q_inv_new,v_inv_new,v_inv_corrected,g_inv,v_f_approx,v_s_f_approx,characteristic_speed,diss_l,diss_r,alpha_loc,sign,bound_number,spec,geom_coeff,geom_source_R,geom_source_Q,phys_source_R,phys_source_Q,phys_source_contact,vel_state,Y_state,S_mom,S_rhoY,g_vel)

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
					r_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
					R_inv_half	= v%pr(dim)%cells(i,j,k)								+ G_half*(p%cells(i,j,k))
					R_inv_old	= v_old(dim,i,j,k)										+ G_half_old*(p_old(i,j,k))
					
					q_inv(1) 	= v_f(dim,dim,i,j,k)									- G_half*p_f(dim,i,j,k)
					q_inv(2) 	= v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	- G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
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
					! balance source.  Optional physical sources use the same effective
					! conservative rates that have already entered the predictor and will
					! enter the corrector; no residual-based g_inv reconstruction is used.
					geom_coeff = 0.0_dp
					if (cabaret_use_geometry_source_shifts .and. dim == 1) then
						geom_coeff = finite_volume_geometry_coefficient(nu, mesh%mesh(1,i,j,k), cell_size(1))
					end if
					geom_source_R = -geom_coeff * v_s%cells(i,j,k) * v%pr(1)%cells(i,j,k)
					geom_source_Q =  geom_coeff * v_s%cells(i,j,k) * v%pr(1)%cells(i,j,k)

					phys_source_R = 0.0_dp
					phys_source_Q = 0.0_dp
					phys_source_contact = 0.0_dp
					g_vel(:) = 0.0_dp
					if (cabaret_use_effective_physical_sources .and. cabaret_use_effective_physical_source_shifts) then
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
					end if

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
						v_inv_corr(1,dim2,dim,i,j,k) = limit_quasi_invariant(v_inv_new(dim2,1), v_inv(dim2,1), v_inv_half(dim2), v_inv(dim2,2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
						v_inv_corr(2,dim2,dim,i,j,k) = limit_quasi_invariant(v_inv_new(dim2,2), v_inv(dim2,1), v_inv_half(dim2), v_inv(dim2,2), &
															  g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
		
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

		deallocate(v_inv, v_inv_new, v_inv_half, v_inv_old)
		deallocate(vel_state, S_mom, g_vel, Y_state, S_rhoY)
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
            
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,loop,G_half_lower,G_half_higher,v_f_approx,v_s_f_approx,characteristic_speed,sign,bound_number)

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
							rho_f_new%cells(dim,i,j,k)			=	density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
											    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								
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
							rho_f_new%cells(dim,i,j,k)			=	density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
							
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
							rho_f_new%cells(dim,i,j,k)	= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
											    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
						
							do dim1 = 1,dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	0.5_dp * (r_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + q_inv_corr(2,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(2,dim1,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end if
							end do	
						end if
					
						if (( characteristic_speed(1) < 0.0_dp ).and.&
							( characteristic_speed(2) < 0.0_dp ).and.&
							( characteristic_speed(3) < 0.0_dp )) then
						
							p_f_new%cells(dim,i,j,k)			= 0.5_dp * (r_inv_corr(1,dim,i,j,k) - q_inv_corr(1,dim,i,j,k)) / G_half_higher
							rho_f_new%cells(dim,i,j,k)			= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
						
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
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
										    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							
							do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
                                
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
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(2,dim,dim,i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), &
										    rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							
                            do dim1 = 1,dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corr(1,dim1,dim,i,j,k)
                                end if
							end do
                        else
							rho_f_new%cells(dim,i,j,k)		= density_from_contact_quasi_invariant(p_f_new%cells(dim,i,j,k), v_inv_corr(1,dim,dim,i,j,k), rho%cells(i,j,k), v_s%cells(i,j,k))
                                
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


	subroutine reconstruct_contact_family_face_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, spec, bound_number
		integer, dimension(3,2) :: cons_inner_loop, cons_utter_loop, loop
		real(dp), dimension(3) :: cell_size
		real(dp), dimension(:,:), allocatable :: Y_inv, Y_inv_corrected, Y_inv_new
		real(dp), dimension(:), allocatable :: Y_inv_half, Y_inv_old
		real(dp), dimension(:), allocatable :: vel_state, Y_state, S_mom, S_rhoY, g_Y
		real(dp) :: g_inv, alpha_loc, g_entropy
		real(dp) :: v_f_approx_lower, v_f_approx_higher
		real(dp) :: diss_l, diss_r
		real(dp), dimension(2) :: S_inv, S_inv_new, S_inv_corrected
		real(dp) :: S_inv_half
		real(dp) :: T_face_lower, T_face_higher

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cons_utter_loop = this%domain%get_local_utter_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

		call this%exchange_face_pressure_density()
		call this%initialize_material_entropy_faces()

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
					s_material_f_new => this%s_material_f_new, &
					bc			=> this%boundary%bc_ptr)
        
		!$omp parallel default(shared)  private(i,j,k,dim,loop,spec,Y_inv,Y_inv_half,Y_inv_new,Y_inv_old,g_inv,alpha_loc,Y_inv_corrected,v_f_approx_lower,v_f_approx_higher,bound_number,diss_r,diss_l,S_inv,S_inv_new,S_inv_corrected,S_inv_half,T_face_lower,T_face_higher,vel_state,Y_state,S_mom,S_rhoY,g_Y,g_entropy) !, &
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

					! Entropy is reconstructed with the same contact-family stencil as Y_k.
					! The final density closure will use this entropy together with the
					! normalized face mass fractions, so the density and composition
					! contacts cannot drift apart numerically.
					T_face_lower  = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(this%p_f(dim,i,j,k), &
						this%rho_f(dim,i,j,k), Y_f(:,dim,i,j,k))
					T_face_higher = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(this%p_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)), &
						this%rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)), Y_f(:,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))
					S_inv(1) = this%thermo%thermo_ptr%mixture_specific_entropy(T_face_lower, this%p_f(dim,i,j,k), Y_f(:,dim,i,j,k), cabaret_entropy_p_ref)
					S_inv(2) = this%thermo%thermo_ptr%mixture_specific_entropy(T_face_higher, this%p_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)), &
						Y_f(:,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)), cabaret_entropy_p_ref)
					S_inv_half = this%thermo%thermo_ptr%mixture_specific_entropy(T%cells(i,j,k), p%cells(i,j,k), Y_inv_half(:), cabaret_entropy_p_ref)
					
					diss_l = 0.0_dp
					diss_r = 0.0_dp

					do spec = 1,species_number
						y_inv_new(spec,1)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_l)*y_inv(spec,2))/(1.0_dp+diss_l)
						y_inv_new(spec,2)	= (2.0_dp*Y_inv_half(spec) - (1.0_dp-diss_r)*y_inv(spec,1))/(1.0_dp+diss_r)
					end do	

					S_inv_new(1) = (2.0_dp*S_inv_half - (1.0_dp-diss_l)*S_inv(2))/(1.0_dp+diss_l)
					S_inv_new(2) = (2.0_dp*S_inv_half - (1.0_dp-diss_r)*S_inv(1))/(1.0_dp+diss_r)

					g_Y(:) = 0.0_dp
					g_entropy = 0.0_dp
					if (cabaret_use_effective_physical_sources .and. cabaret_use_effective_physical_source_shifts) then
						do spec = 1, species_number
							Y_state(spec) = Y%pr(spec)%cells(i,j,k)
							S_rhoY(spec) = rhoY_src_old(spec,i,j,k)
						end do
						do spec = 1, dimensions
							vel_state(spec) = v%pr(spec)%cells(i,j,k)
							S_mom(spec) = mom_src_old(spec,i,j,k)
						end do
						call this%physical_source_contact_shifts(rho%cells(i,j,k), p%cells(i,j,k), v_s%cells(i,j,k), &
							E_f%cells(i,j,k), vel_state, Y_state, rho_src_old(i,j,k), S_mom, &
							rhoE_src_old(i,j,k), S_rhoY, g_Y, g_entropy)
					end if

					g_inv = g_entropy
					alpha_loc = 0.0_dp
					S_inv_corrected(1) = limit_quasi_invariant(S_inv_new(1), S_inv(1), S_inv_half, S_inv(2), &
						g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					S_inv_corrected(2) = limit_quasi_invariant(S_inv_new(2), S_inv(1), S_inv_half, S_inv(2), &
						g_inv*this%time_step, alpha_loc, -huge(1.0_dp), huge(1.0_dp))
					
					do spec = 1,species_number
						! Geometry changes rho and rho*Y_k by the same factor; physical
						! effective sources use dY_k/dt=(S_rhoY-Y_k*S_rho)/rho.
						g_inv = g_Y(spec)
						alpha_loc = 0.0_dp
						y_inv_corrected(spec,1) = limit_quasi_invariant(y_inv_new(spec,1), y_inv(spec,1), Y_inv_half(spec), y_inv(spec,2), &
													  g_inv*this%time_step, alpha_loc, 0.0_dp, 1.0_dp)
						y_inv_corrected(spec,2) = limit_quasi_invariant(y_inv_new(spec,2), y_inv(spec,1), Y_inv_half(spec), y_inv(spec,2), &
													  g_inv*this%time_step, alpha_loc, 0.0_dp, 1.0_dp)
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
								! Consistency fix: upwind Y from the same reconstructed face velocity as the gas-dynamic state.
								v_f_approx_lower		= v_f_new%pr(dim)%cells(dim,i,j,k)
								if (v_f_approx_lower < 0.0_dp) then
                                	Y_f_new%pr(spec)%cells(dim,i,j,k) =  (y_inv_corrected(spec,1))
									if (spec == 1) s_material_f_new(dim,i,j,k) = S_inv_corrected(1)
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
									if (spec == 1) s_material_f_new(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = S_inv_corrected(2)
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


	subroutine correct_conservative_full_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dimensions, species_number, nu
		integer :: i, j, k, dim, dim1, spec
		integer, dimension(3,2) :: cons_inner_loop
		real(dp), dimension(3) :: cell_size
		character(len=20) :: coordinate_system
		real(dp) :: spec_summ, mean_higher, mean_lower, r
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
            
		!$omp parallel default(shared)  private(i,j,k,dim,dim1,spec,spec_summ,mean_higher,mean_lower,r,rhoY_new)
		allocate(rhoY_new(species_number))
         
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

				if (cabaret_use_effective_physical_sources) then
					rho%cells(i,j,k) = rho_old(i,j,k) + 0.5_dp*this%time_step * &
						(rho_rhs_old(i,j,k) + rho%cells(i,j,k)) + this%time_step*rho_src_old(i,j,k)
				else
					rho%cells(i,j,k)	=	rho_old(i,j,k) + 0.5_dp*this%time_step * &
						(rho%cells(i,j,k) + rho_src_old(i,j,k))
				end if
	   
				! Species are advanced as conservative partial densities rho*Y_k.
				! Mass fractions are recovered after positivity correction of rho*Y_k.
				spec_summ = 0.0_dp
				do spec = 1,species_number
					rhoY_new(spec) = 0.0_dp
					do	dim = 1,dimensions
						mean_higher	=  rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  * y_f_new%pr(spec)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))
						mean_lower	=  rho_f_new%cells(dim,i,j,k)    * y_f_new%pr(spec)%cells(dim,i,j,k)   *v_f_new%pr(dim)%cells(dim,i,j,k)
						rhoY_new(spec)	=  rhoY_new(spec) - (mean_higher - mean_lower ) /cell_size(1)
						if(dim == 1) then
							rhoY_new(spec)	=  rhoY_new(spec) - 2.0_dp * (nu - 1)/((r + 0.5_dp*cell_size(1))**(nu - 1) + (r - 0.5_dp*cell_size(1))**(nu - 1))		&
																	* 0.5_dp * (mean_higher + mean_lower)														&
																	* ((r + 0.5_dp*cell_size(1))**(nu - 1) - (r - 0.5_dp*cell_size(1))**(nu - 1))											&
																	/ ((r + 0.5_dp*cell_size(1)) - (r - 0.5_dp*cell_size(1)))
						end if
					end do
					if (cabaret_use_effective_physical_sources) then
						rhoY_new(spec) = rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * &
							(rhoY_rhs_old(spec,i,j,k) + rhoY_new(spec)) + this%time_step*rhoY_src_old(spec,i,j,k)
					else
						rhoY_new(spec) = rho_old(i,j,k) * Y_old(spec,i,j,k) + 0.5_dp * this%time_step * &
							(rhoY_new(spec) + rhoY_src_old(spec,i,j,k))
					end if
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
					
					if (cabaret_use_effective_physical_sources) then
						v%pr(dim)%cells(i,j,k) = rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dp*this%time_step * &
							(mom_rhs_old(dim,i,j,k) + v%pr(dim)%cells(i,j,k)) + this%time_step*mom_src_old(dim,i,j,k)
					else
						v%pr(dim)%cells(i,j,k)	= rho_old(i,j,k)*v_old(dim,i,j,k) +  0.5_dp * this%time_step * &
							(v%pr(dim)%cells(i,j,k) + mom_src_old(dim,i,j,k))
					end if
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
				if (cabaret_use_effective_physical_sources) then
					E_f%cells(i,j,k) = rho_old(i,j,k) * E_f_old(i,j,k) + 0.5_dp*this%time_step * &
						(rhoE_rhs_old(i,j,k) + E_f%cells(i,j,k)) + this%time_step*rhoE_src_old(i,j,k)
				else
					E_f%cells(i,j,k) = rho_old(i,j,k) * E_f_old(i,j,k) + 0.5_dp * this%time_step * &
						(E_f%cells(i,j,k) + rhoE_src_old(i,j,k))
				end if
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


	subroutine finalize_gas_dynamics_step(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: particles_phase_counter

        associate(	v   => this%v%v_ptr)
			call this%mpi_support%exchange_conservative_vector_field(v)	
        end associate

		call this%apply_boundary_conditions_main()

		if(this%additional_particles_phases_number /= 0 .and. &
			(.not. (cabaret_use_effective_physical_sources .and. cabaret_include_particles_in_gas_step))) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
    			call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_main(this%time)
			end do
		end if
	end subroutine finalize_gas_dynamics_step


	subroutine solve_split_physics(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: particles_phase_counter

        call cabaret_heattransfer_timer%tic()
		if (this%heat_trans_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_heat_transfer_in_gas_step))) &
			call this%heat_trans_solver%solve_heat_transfer(this%time_step)
        call cabaret_heattransfer_timer%toc(new_iter=.true.)
        
        call cabaret_diffusion_timer%tic()
		if (this%diffusion_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_diffusion_in_gas_step))) &
			call this%diff_solver%solve_diffusion(this%time_step)
        call cabaret_diffusion_timer%toc(new_iter=.true.)
        
        call cabaret_viscosity_timer%tic()
		if (this%viscosity_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_viscosity_in_gas_step))) &
			call this%viscosity_solver%solve_viscosity(this%time_step)
        call cabaret_viscosity_timer%toc(new_iter=.true.)
        
        call cabaret_chemistry_timer%tic()
		if (this%reactive_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_chemistry_in_gas_step))) &
			call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
        call cabaret_chemistry_timer%toc(new_iter=.true.)
        

		if(this%additional_particles_phases_number /= 0 .and. &
			(.not. (cabaret_use_effective_physical_sources .and. cabaret_include_particles_in_gas_step))) then
			do particles_phase_counter = 1, this%additional_particles_phases_number
				call this%particles_solver(particles_phase_counter)%particles_solve(this%time_step)				!# Lagrangian particles solver
!				call this%particles_solver(particles_phase_counter)%particles_euler_step_v_E(this%time_step)	!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%apply_boundary_conditions_interm_v_d()		!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%particles_lagrange_step(this%time_step)		!# Continuum particles solver
!				call this%particles_solver(particles_phase_counter)%particles_final_step(this%time_step)		!# Continuum particles solver		
			end do		
		end if 	        
	end subroutine solve_split_physics


	subroutine apply_split_sources(this)
		class(cabaret_solver), intent(inout) :: this

		real(dp) :: energy_output_time = 2.0e-07_dp
		real(dp) :: energy_output_radii = 4.0e-04_dp
		real(dp) :: energy_source = 1.9e+04_dp
		real(dp), save :: energy_output_rho = 0.0_dp
		real(dp), save :: energy_output = 0.0_dp
		integer, save :: energy_output_flag = 0
		integer :: dimensions, species_number
		integer :: i, j, k, dim, spec, particles_phase_counter
		integer, dimension(3,2) :: cons_inner_loop
		real(dp), dimension(3) :: cell_size
		real(dp) :: spec_summ

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()
		cell_size       = this%mesh%mesh_ptr%get_cell_edges_length()

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
		
		!$omp parallel default(shared)  private(i,j,k,dim,spec,spec_summ)

		!$omp do collapse(3) schedule(guided)		
		do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
			
			if(bc%bc_markers(i,j,k) == 0) then	
				
				if (this%reactive_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_chemistry_in_gas_step)))	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_chem%cells(i,j,k) * this%time_step
					do spec = 1,species_number
						Y_prod(spec,i,j,k)	= Y_prod(spec,i,j,k)	+ Y_prod_chem%pr(spec)%cells(i,j,k) * this%time_step
					end do		
                end if
                
				if (this%heat_trans_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_heat_transfer_in_gas_step)))	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_heat%cells(i,j,k) * this%time_step
                end if
                
 				if (this%diffusion_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_diffusion_in_gas_step)))	then
					E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_diff%cells(i,j,k) * this%time_step
					do spec = 1, species_number
						Y_prod(spec,i,j,k)	= Y_prod(spec,i,j,k)	+ Y_prod_diff%pr(spec)%cells(i,j,k) * this%time_step
					end do
				end if

				if (this%viscosity_flag .and. (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_viscosity_in_gas_step)))	then
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
                
                if (this%additional_particles_phases_number /= 0 .and. &
                    (.not. (cabaret_use_effective_physical_sources .and. cabaret_include_particles_in_gas_step))) then
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
	end subroutine apply_split_sources


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

	subroutine exchange_face_pressure_density(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
	end subroutine exchange_face_pressure_density

	subroutine exchange_face_primitive_state(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_primitive_state

	subroutine exchange_face_thermodynamic_state(this)
		class(cabaret_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%e_i_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%E_f_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%v_s_f_new%s_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_thermodynamic_state

	subroutine zero_new_flow_state(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, spec, dimensions, species_number

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number

		this%p_f_new%s_ptr%cells   = 0.0_dp
		this%rho_f_new%s_ptr%cells = 0.0_dp
		if (allocated(this%s_material_f_new)) this%s_material_f_new = 0.0_dp

		do dim = 1, dimensions
			this%v_f_new%v_ptr%pr(dim)%cells = 0.0_dp
		end do

		do spec = 1, species_number
			this%Y_f_new%v_ptr%pr(spec)%cells = 0.0_dp
		end do
	end subroutine zero_new_flow_state


	subroutine initialize_material_entropy_faces(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, dim1, i, j, k, dimensions
		integer, dimension(3,2) :: flow_inner_loop, loop
		real(dp) :: T_face

		if (.not. cabaret_use_entropy_locked_contact_density) return

		dimensions      = this%domain%get_domain_dimensions()
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()

		associate( p_f => this%p_f, rho_f => this%rho_f, Y_f => this%Y_f, s_f => this%s_material_f_new )
		!$omp parallel default(shared) private(dim,dim1,i,j,k,loop,T_face)
		do dim = 1, dimensions
			loop = flow_inner_loop
			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

			!$omp do collapse(3) schedule(static)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				T_face = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_f(dim,i,j,k), rho_f(dim,i,j,k), Y_f(:,dim,i,j,k))
				s_f(dim,i,j,k) = this%thermo%thermo_ptr%mixture_specific_entropy(T_face, p_f(dim,i,j,k), Y_f(:,dim,i,j,k), cabaret_entropy_p_ref)
			end do
			end do
			end do
			!$omp end do nowait
		end do
		!$omp end parallel
		end associate
	end subroutine initialize_material_entropy_faces


	subroutine enforce_material_contact_density_from_entropy(this)
		class(cabaret_solver), intent(inout) :: this

		integer :: dim, dim1, i, j, k, dimensions, spec, species_number
		integer, dimension(3,2) :: flow_inner_loop, loop
		integer :: il, jl, kl
		integer :: temperature_guard_count, molar_volume_guard_count
		real(dp) :: m_left, m_right, molar_ratio, T_guess, T_entropy, T_used, T_guard
		real(dp) :: T_left, T_right, T_ref, m_face
		real(dp) :: molar_denom_left, molar_denom_right
		real(dp) :: psi_left, psi_right, psi_min, psi_max, psi_span
		real(dp) :: psi_clip_tol, psi_trigger_tol, psi_entropy, psi_used
		real(dp) :: psi_lower_clip, psi_upper_clip, psi_lower_trigger, psi_upper_trigger
		real(dp) :: T_consequence_limit
		logical :: apply_entropy_lock, apply_temperature_guard, apply_molar_volume_guard
		logical :: strong_guard_contact, temperature_consequence
		real(dp), dimension(:), allocatable :: Y_face

		if (.not. cabaret_use_entropy_locked_contact_density) return

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number
		flow_inner_loop = this%domain%get_local_inner_faces_bounds()
		temperature_guard_count   = 0
		molar_volume_guard_count = 0

		associate( rho_f_new => this%rho_f_new%s_ptr, &
				   p_f_new   => this%p_f_new%s_ptr, &
				   Y_f_new   => this%Y_f_new%v_ptr, &
				   Y         => this%Y%v_ptr, &
				   rho       => this%rho%s_ptr, &
				   T         => this%T%s_ptr, &
				   s_f       => this%s_material_f_new, &
				   bc        => this%boundary%bc_ptr )

		!$omp parallel default(shared) private(dim,dim1,i,j,k,loop,il,jl,kl,spec,Y_face, &
		!$omp& m_left,m_right,molar_ratio,T_guess,T_entropy,T_used,T_guard,T_left,T_right,T_ref, &
		!$omp& m_face,molar_denom_left,molar_denom_right,psi_left,psi_right,psi_min,psi_max, &
		!$omp& psi_span,psi_clip_tol,psi_trigger_tol,psi_entropy,psi_used,psi_lower_clip, &
		!$omp& psi_upper_clip,psi_lower_trigger,psi_upper_trigger,T_consequence_limit, &
		!$omp& apply_entropy_lock,apply_temperature_guard,apply_molar_volume_guard, &
		!$omp& strong_guard_contact,temperature_consequence) &
		!$omp& reduction(+:temperature_guard_count,molar_volume_guard_count)
		allocate(Y_face(species_number))
		do dim = 1, dimensions
			loop = flow_inner_loop
			do dim1 = 1, dimensions
				loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
			end do

			!$omp do collapse(3) schedule(guided)
			do k = loop(3,1), loop(3,2)
			do j = loop(2,1), loop(2,2)
			do i = loop(1,1), loop(1,2)
				il = i - I_m(dim,1)
				jl = j - I_m(dim,2)
				kl = k - I_m(dim,3)

				if ((bc%bc_markers(i,j,k) == 0) .and. (bc%bc_markers(il,jl,kl) == 0)) then
					molar_denom_left  = 0.0_dp
					molar_denom_right = 0.0_dp
					do spec = 1, species_number
						molar_denom_left  = molar_denom_left  + Y%pr(spec)%cells(il,jl,kl)/this%thermo%thermo_ptr%molar_masses(spec)
						molar_denom_right = molar_denom_right + Y%pr(spec)%cells(i,j,k)/this%thermo%thermo_ptr%molar_masses(spec)
					end do
					m_left  = 1.0_dp/molar_denom_left
					m_right = 1.0_dp/molar_denom_right
					molar_ratio = max(m_left,m_right)/min(m_left,m_right)

					apply_entropy_lock = cabaret_entropy_lock_all_interior_faces .or. &
						(molar_ratio > cabaret_entropy_lock_molar_mass_ratio)

					if (apply_entropy_lock) then
						do spec = 1, species_number
							Y_face(spec) = Y_f_new%pr(spec)%cells(dim,i,j,k)
						end do

						T_guess   = this%thermo%thermo_ptr%temperature_from_pressure_density_Y(p_f_new%cells(dim,i,j,k), rho_f_new%cells(dim,i,j,k), Y_face)
						T_entropy = this%thermo%thermo_ptr%temperature_from_entropy_pressure_Y(s_f(dim,i,j,k), p_f_new%cells(dim,i,j,k), Y_face, T_guess, &
			cabaret_entropy_temperature_bracket_low, cabaret_entropy_temperature_bracket_high, cabaret_entropy_p_ref)

						T_used = T_entropy
						apply_temperature_guard = .false.
						if (cabaret_use_entropy_lock_temperature_guard) then
							T_left  = T%cells(il,jl,kl)
							T_right = T%cells(i,j,k)
							T_ref   = max(T_guess, T_left, T_right)
							T_guard = T_ref*cabaret_entropy_lock_temperature_guard_ratio + &
								cabaret_entropy_lock_temperature_guard_abs
							apply_temperature_guard = (T_entropy > T_guard)
							if (apply_temperature_guard) then
								T_used = T_guard
								temperature_guard_count = temperature_guard_count + 1
							end if
						end if

						m_face = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_face)

						! Entropy locking gives psi = M/rho = R_u*T/p.  The remaining
						! one-cell spike is a non-monotone psi outlier: composition is
						! already close to the heavy-gas side while density is still too
						! close to the light-gas/contact-transition side.  Clip psi, not
						! p_f, Y_f, s_f or velocities.
						psi_used = r_gase_J*T_used/ &
							p_f_new%cells(dim,i,j,k)

						apply_molar_volume_guard = .false.
						if (cabaret_use_entropy_lock_molar_volume_guard) then
							strong_guard_contact = cabaret_entropy_lock_all_interior_faces .or. &
								(molar_ratio > cabaret_entropy_lock_molar_volume_guard_molar_mass_ratio)

							if (strong_guard_contact) then
								psi_left  = m_left /rho%cells(il,jl,kl)
								psi_right = m_right/rho%cells(i,j,k)
								psi_min   = min(psi_left,psi_right)
								psi_max   = max(max(psi_left,psi_right), psi_min)
								psi_span  = psi_max - psi_min
								psi_clip_tol = cabaret_entropy_lock_molar_volume_guard_clip_abs_tol + &
									cabaret_entropy_lock_molar_volume_guard_clip_rel_tol*psi_span
								psi_trigger_tol = cabaret_entropy_lock_molar_volume_guard_trigger_abs_tol + &
									cabaret_entropy_lock_molar_volume_guard_trigger_rel_tol*psi_span
								psi_lower_clip    = psi_min - psi_clip_tol
								psi_upper_clip    = max(psi_lower_clip, psi_max + psi_clip_tol)
								psi_lower_trigger = psi_min - psi_trigger_tol
								psi_upper_trigger = max(psi_lower_trigger, psi_max + psi_trigger_tol)
								psi_entropy = psi_used

								T_left  = T%cells(il,jl,kl)
								T_right = T%cells(i,j,k)
								T_ref   = max(T_guess, T_left, T_right)
								T_consequence_limit = T_ref*(1.0_dp + cabaret_entropy_lock_molar_volume_guard_T_ratio) + &
									cabaret_entropy_lock_molar_volume_guard_T_abs
								temperature_consequence = (T_entropy > T_consequence_limit)

								if (psi_entropy <= 0.0_dp) then
									write(*,*) 'CABARET strict molar-volume guard: invalid psi at ', dim, i, j, k
									write(*,*) '  psi_entropy, T_used, p_f = ', psi_entropy, T_used, p_f_new%cells(dim,i,j,k)
									error stop 'CABARET strict molar-volume guard: invalid psi'
								end if

								if (temperature_consequence .and. (psi_entropy > psi_upper_trigger)) then
									psi_used = psi_upper_clip
									apply_molar_volume_guard = .true.
								else if ((.not. cabaret_entropy_lock_molar_volume_guard_upper_only) .and. &
									temperature_consequence .and. (psi_entropy < psi_lower_trigger)) then
									psi_used = psi_lower_clip
									apply_molar_volume_guard = .true.
								end if

								if (apply_molar_volume_guard) then
									molar_volume_guard_count = molar_volume_guard_count + 1
								end if
							end if
						end if

						rho_f_new%cells(dim,i,j,k) = m_face/psi_used
					end if
				end if
			end do
			end do
			end do
			!$omp end do nowait
		end do
		deallocate(Y_face)
		!$omp end parallel
		end associate

		if (cabaret_entropy_lock_temperature_guard_print_statistics .and. temperature_guard_count > 0) then
			print *, 'CABARET entropy-lock temperature guard clipped faces:', temperature_guard_count
		end if
		if (cabaret_entropy_lock_molar_volume_guard_print_statistics .and. molar_volume_guard_count > 0) then
			print *, 'CABARET entropy-lock molar-volume guard clipped faces:', molar_volume_guard_count
		end if

		! rho_f_new has changed after the regular face exchange in finish_face_reconstruction().
		call this%exchange_face_pressure_density()
	end subroutine enforce_material_contact_density_from_entropy

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
								this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)), 0.0_dp)
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


	real(dp) function cell_length(this, dim) result(dx)
		class(cabaret_solver), intent(in) :: this
		integer, intent(in) :: dim
		real(dp), dimension(3) :: cell_size

		cell_size = this%mesh%mesh_ptr%get_cell_edges_length()
		if (dim >= 1 .and. dim <= size(cell_size) .and. cell_size(dim) > 0.0_dp) then
			dx = cell_size(dim)
		else
			dx = cell_size(1)
		end if
	end function cell_length

	logical function state_is_finite(this, value) result(ok)
		class(cabaret_solver), intent(in) :: this
		real(dp), intent(in) :: value

		ok = (abs(value) < huge(1.0_dp))
	end function state_is_finite

	subroutine check_conservative_state(this, stage_name)
		class(cabaret_solver), intent(in) :: this
		character(len=*), intent(in) :: stage_name

		integer :: i, j, k, spec, species_number
		integer, dimension(3,2) :: cons_inner_loop
		real(dp) :: y_sum, kinetic, e_internal

		if (.not. cabaret_debug_checks) return

		species_number  = this%chem%chem_ptr%species_number
		cons_inner_loop = this%domain%get_local_inner_cells_bounds()

		do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
		do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
		do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
			if (this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then
				if (this%rho%s_ptr%cells(i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET state check: non-positive rho at ', trim(stage_name), i, j, k, this%rho%s_ptr%cells(i,j,k)
				end if
				if (this%p%s_ptr%cells(i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET state check: non-positive p at ', trim(stage_name), i, j, k, this%p%s_ptr%cells(i,j,k)
				end if
				if (this%T%s_ptr%cells(i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET state check: non-positive T at ', trim(stage_name), i, j, k, this%T%s_ptr%cells(i,j,k)
				end if

				y_sum = 0.0_dp
				do spec = 1, species_number
					y_sum = y_sum + this%Y%v_ptr%pr(spec)%cells(i,j,k)
				end do
				if (abs(y_sum - 1.0_dp) > 1.0e-8_dp) then
					write(*,*) 'CABARET state check: sum(Y) /= 1 at ', trim(stage_name), i, j, k, y_sum
				end if

				kinetic = 0.0_dp
				! The loop over velocity components is intentionally written explicitly through size guards.
				if (this%domain%get_domain_dimensions() >= 1) kinetic = kinetic + 0.5_dp*this%v%v_ptr%pr(1)%cells(i,j,k)**2
				if (this%domain%get_domain_dimensions() >= 2) kinetic = kinetic + 0.5_dp*this%v%v_ptr%pr(2)%cells(i,j,k)**2
				if (this%domain%get_domain_dimensions() >= 3) kinetic = kinetic + 0.5_dp*this%v%v_ptr%pr(3)%cells(i,j,k)**2
				e_internal = this%E_f%s_ptr%cells(i,j,k) - kinetic
				if (e_internal <= 0.0_dp) then
					write(*,*) 'CABARET state check: non-positive internal energy at ', trim(stage_name), i, j, k, e_internal
				end if
			end if
		end do
		end do
		end do
	end subroutine check_conservative_state

	subroutine check_face_state(this, stage_name)
		class(cabaret_solver), intent(in) :: this
		character(len=*), intent(in) :: stage_name

		integer :: i, j, k, dim, dim1, spec, dimensions, species_number
		integer, dimension(3,2) :: flow_inner_loop, loop
		real(dp) :: y_sum, kinetic, e_internal

		if (.not. cabaret_debug_checks) return

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
				if (this%rho_f_new%s_ptr%cells(dim,i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET face check: non-positive rho_f at ', trim(stage_name), dim, i, j, k, this%rho_f_new%s_ptr%cells(dim,i,j,k)
				end if
				if (this%p_f_new%s_ptr%cells(dim,i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET face check: non-positive p_f at ', trim(stage_name), dim, i, j, k, this%p_f_new%s_ptr%cells(dim,i,j,k)
				end if
				if (this%T_f_new%s_ptr%cells(dim,i,j,k) <= 0.0_dp) then
					write(*,*) 'CABARET face check: non-positive T_f at ', trim(stage_name), dim, i, j, k, this%T_f_new%s_ptr%cells(dim,i,j,k)
				end if

				y_sum = 0.0_dp
				do spec = 1, species_number
					y_sum = y_sum + this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
				end do
				if (abs(y_sum - 1.0_dp) > 1.0e-8_dp) then
					write(*,*) 'CABARET face check: sum(Y_f) /= 1 at ', trim(stage_name), dim, i, j, k, y_sum
				end if

				kinetic = 0.0_dp
				do dim1 = 1, dimensions
					kinetic = kinetic + 0.5_dp*this%v_f_new%v_ptr%pr(dim1)%cells(dim,i,j,k)**2
				end do
				e_internal = this%E_f_f_new%s_ptr%cells(dim,i,j,k) - kinetic
				if (e_internal <= 0.0_dp) then
					write(*,*) 'CABARET face check: non-positive internal energy at ', trim(stage_name), dim, i, j, k, e_internal
				end if
			end do
			end do
			end do
		end do
	end subroutine check_face_state

	subroutine apply_boundary_conditions_flow(this, dim,i,j,k, characteristic_speed, q_inv_corrected, r_inv_corrected, v_inv_corrected, G_half)

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
					mol_mix_conc	=> this%mol_mix_conc%s_ptr		, &
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
										mol_mix_conc%cells(ghost_i,ghost_j,ghost_k)	= mol_mix_conc%cells(i,j,k)
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
										mol_mix_conc%cells(ghost_i,ghost_j,ghost_k)	= mol_mix_conc%cells(i,j,k)
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
											mol_mix_conc%cells(ghost_i,ghost_j,ghost_k) = mol_mix_from_farfield
										else
											mol_mix_conc%cells(ghost_i,ghost_j,ghost_k) = mol_mix_conc%cells(i,j,k)
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
											mol_mix_conc%cells(ghost_i,ghost_j,ghost_k)	= mol_mix_conc%cells(i,j,k)
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
												mol_mix_conc%cells(ghost_i,ghost_j,ghost_k) = mol_mix_from_farfield
											else
												mol_mix_conc%cells(ghost_i,ghost_j,ghost_k) = mol_mix_conc%cells(i,j,k)
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
	
	subroutine calculate_time_step(this)

#ifdef mpi
	use MPI
#endif

		class(cabaret_solver)	,intent(inout)	:: this
		
		real(dp)	:: local_min_dt, global_min_dt, delta_t_interm, velocity_value
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
				velocity_value = 0.0_dp
				directional_wave_sum = 0.0_dp
				do dim = 1,dimensions
					velocity_value = velocity_value + v%pr(dim)%cells(i,j,k)*v%pr(dim)%cells(i,j,k)
					if (cell_size(dim) > 0.0_dp) then
						directional_wave_sum = directional_wave_sum + (abs(v%pr(dim)%cells(i,j,k)) + v_s%cells(i,j,k)) / cell_size(dim)
					end if
				end do

				if (directional_wave_sum > 0.0_dp) then
					delta_t_interm = 1.0_dp / directional_wave_sum
				else
					delta_t_interm = minval(cell_size, mask = cell_size > 0.0_dp) / sqrt(velocity_value) + v_s%cells(i,j,k)
				end if

				local_min_dt = min(local_min_dt, delta_t_interm)
			end if
		end do
		end do
		end do
	
#ifdef mpi					
		call mpi_allreduce(local_min_dt, global_min_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpi_communicator, error)
#else
		global_min_dt = local_min_dt
#endif

		if (global_min_dt < huge(1.0_dp)) then
			this%time_step = this%courant_fraction * global_min_dt
		else
			! Keep the previous value if all local cells are inactive.
			this%time_step = this%time_step
		end if

		end associate
			
	end subroutine calculate_time_step

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

	!subroutine check_symmetry(this)
 !   
	!	class(cabaret_solver)	,intent(inout)	:: this
	!	
	!	real(dp)	:: delta_t_interm, time_step(1), velocity_value
	!	real(dp)	,dimension(:)	,allocatable	,save	:: time_step_array
 !
	!	integer						:: dimensions, species_number
	!	integer						:: processor_rank, processor_number, mpi_communicator
 !
	!	integer		,dimension(3,2)	:: cons_inner_loop, cons_utter_loop, flow_inner_loop
 !       integer		,dimension(3,2)	:: loop
	!	real(dp)	,dimension(3)	:: cell_size
	!	integer	:: sign
	!	integer :: i,j,k,dim,dim1,error, spec
 !
	!	associate(	rho			=> this%rho%s_ptr		, &
	!				p			=> this%p%s_ptr			, &
	!				E_f			=> this%E_f%s_ptr		, &
	!				v			=> this%v%v_ptr			, &
	!				Y			=> this%Y%v_ptr			, &
	!				
	!				v_f			=> this%v_f		, &
	!				rho_f		=> this%rho_f	, &
	!				E_f_f		=> this%E_f_f	, &
	!				e_i_f		=> this%e_i_f	, &
	!				p_f			=> this%p_f		, &
	!				v_s_f		=> this%v_s_f	, & 
	!				Y_f			=> this%Y_f		, &
	!				
	!				p_f_new		=> this%p_f_new%s_ptr		, &	
	!				rho_f_new	=> this%rho_f_new%s_ptr		, &
	!				E_f_f_new	=> this%E_f_f_new%s_ptr		, &
	!				Y_f_new		=> this%Y_f_new%v_ptr		, &
	!				v_f_new		=> this%v_f_new%v_ptr		, &
 !
	!				bc				=> this%boundary%bc_ptr		, &
	!				mesh			=> this%mesh%mesh_ptr)
	!	
	!	dimensions			= this%domain%get_domain_dimensions()
	!	species_number      = this%chem%chem_ptr%species_number
	!	cons_inner_loop		= this%domain%get_local_inner_cells_bounds()
 !       cons_utter_loop		= this%domain%get_local_utter_cells_bounds()
 !       flow_inner_loop		= this%domain%get_local_inner_faces_bounds()
 !
	!	cell_size			= mesh%get_cell_edges_length()					
	!				
 !       !# central symmetry
 !       
	!	do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
	!	do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
	!	do i = cons_inner_loop(1,1),cons_inner_loop(1,2)
	!		if(bc%bc_markers(i,j,k) == 0) then
	!			if (rho%cells(i,j,k) - rho%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
 !                   print *, "density asymmetry", i,j,k, rho%cells(i,j,k) - rho%cells(cons_inner_loop(1,2)-i,cons_inner_loop(2,2)-j,cons_inner_loop(3,2)-k)
 !                   pause
 !               end if
 !               if (E_f%cells(i,j,k) - E_f%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
 !                   print *, "energy asymmetry", i,j,k, E_f%cells(i,j,k) - E_f%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1)
 !                   pause
 !               end if
 !               do dim = 1, dimensions
	!				if (abs(v%pr(dim)%cells(i,j,k)) - abs(v%pr(dim)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1)) /= 0 ) then
	!					print *, "velocity ",dim, " asymmetry", i,j,k, abs(v%pr(dim)%cells(i,j,k)) - abs(v%pr(dim)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1))
 !                       pause
 !                   end if 
	!			end do
 !               do spec = 1, species_number
	!				if (Y%pr(spec)%cells(i,j,k) - Y%pr(spec)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) /= 0 ) then
	!					print *, "Concentration ",spec, " asymmetry", i,j,k, Y%pr(spec)%cells(i,j,k) - Y%pr(spec)%cells(cons_inner_loop(1,2)-i+1,cons_inner_loop(2,2)-j+1,cons_inner_loop(3,2)-k+1) 
 !                       pause
 !                   end if 
	!			end do            
	!		end if
	!	end do
	!	end do
 !       end do
 !
	!	do dim = 1,dimensions
 !
	!		! Avoid looping in transverse direction in ghost cells
 !
	!		loop(3,1) = cons_inner_loop(3,1)*I_m(dim,3) + cons_inner_loop(3,1)*(1 - I_m(dim,3))
	!		loop(3,2) = cons_utter_loop(3,2)*I_m(dim,3) + cons_inner_loop(3,2)*(1 - I_m(dim,3))
 !
	!		loop(2,1) = cons_inner_loop(2,1)*I_m(dim,2) + cons_inner_loop(2,1)*(1 - I_m(dim,2))
	!		loop(2,2) = cons_utter_loop(2,2)*I_m(dim,2) + cons_inner_loop(2,2)*(1 - I_m(dim,2))	
 !
	!		loop(1,1) = cons_inner_loop(1,1)*I_m(dim,1) + cons_inner_loop(1,1)*(1 - I_m(dim,1))
	!		loop(1,2) = cons_utter_loop(1,2)*I_m(dim,1) + cons_inner_loop(1,2)*(1 - I_m(dim,1))						
 !
	!		!$omp do collapse(3) schedule(guided)	
	!		do k = loop(3,1),loop(3,2)
	!		do j = loop(2,1),loop(2,2)
	!		do i = loop(1,1),loop(1,2)
 !                   
	!			if (rho_f_new%cells(dim,i,j,k) - rho_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
	!				print *, "flow density asymmetry" 
 !                   print *, i,j,k, loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1
 !                   print *, rho_f_new%cells(dim,i,j,k) - rho_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
	!				pause
	!			end if
	!			if (p_f_new%cells(dim,i,j,k) - p_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
	!				print *, "flow pressure asymmetry", i,j,k, p_f_new%cells(dim,i,j,k) - p_f_new%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
	!				pause
	!			end if
	!			do dim1 = 1, dimensions
	!				if (abs(v_f_new%pr(dim)%cells(dim,i,j,k)) - abs(v_f_new%pr(dim)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)) /= 0 ) then
	!					print *, "Flow velocity ",dim, " asymmetry", i,j,k, abs(v_f_new%pr(dim)%cells(dim,i,j,k)) - abs(v_f_new%pr(dim)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1))
	!					pause
	!				end if 
	!			end do
	!			do spec = 1, species_number
	!				if (Y_f_new%pr(spec)%cells(dim,i,j,k) - Y_f_new%pr(spec)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1) /= 0 ) then
	!					print *, "Flow concentration ",spec, " asymmetry", i,j,k, Y_f_new%pr(spec)%cells(dim,i,j,k) - Y_f_new%pr(spec)%cells(dim,loop(1,2)-i+1,loop(2,2)-j+1,loop(3,2)-k+1)
	!					pause
	!				end if 
 !               end do
 !
	!		end do
	!		end do
	!		end do        
 !       
 !       end do
 !       
	!	end associate    
 !   
 !   
 !   
 !   end subroutine
    
  !  subroutine if_stabilized(this,time,stabilized)
		!class(cabaret_solver)	,intent(inout)	:: this
		!real(dp)			,intent(in)			:: time    
  !
		!logical				,intent(out)	:: stabilized
		!
  !      logical						:: boundary 
		!real(dp)	,dimension(3)	:: cell_size		
		!
		!real(dp)					:: max_val, left_val, right_val, flame_velocity, flame_surface_length, surface_factor
		!real(dp)					:: a, b 
		!real(dp)					:: time_diff, time_delay, time_stabilization
		!real(dp), save			:: previous_flame_location = 0.0_dp, current_flame_location = 0.0_dp, farfield_velocity = 0.0_dp
		!real(dp), save			:: previous_time = 0.0_dp, current_time = 0.0_dp
  !      real(dp), save			:: av_flame_velocity = 0.0_dp, previous_av_flame_velocity = 0.0_dp
		!real(dp), dimension(20),	save	:: flame_velocity_array = 0.0_dp
		!integer		,save			:: correction = 0, counter = 0
		!integer						:: flame_front_index
		!character(len=200)			:: file_name
		!
		!
		!integer	:: dimensions, species_number
		!integer	,dimension(3,2)	:: cons_inner_loop
  !
		!real(dp)				:: tip_coord, side_coord_x, side_coord_y, T_flame, lp_dist, x_f, min_dist, min_y, max_x
		!integer					:: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound, lp_start, lp_number2 
		!
		!integer :: CO_index, H2O2_index, HO2_index
		!integer	:: bound_number,sign
		!integer :: i,j,k,plus,dim,dim1,spec, lp_index,lp_index2,lp_index3
		!
		!
		!character(len=20)		:: boundary_type_name
		!
		!dimensions		= this%domain%get_domain_dimensions()
		!species_number	= this%chem%chem_ptr%species_number
		!
		!cons_inner_loop	= this%domain%get_local_inner_cells_bounds()
		!		
		!cell_size		= this%mesh%mesh_ptr%get_cell_edges_length()
  !
		!HO2_index		= this%chem%chem_ptr%get_chemical_specie_index('HO2')
		!
		!stabilized = .false.
		!
		!associate (	v				=> this%v%v_ptr				, &
		!			T				=> this%T%s_ptr				, &
		!			Y				=> this%Y%v_ptr				, &
		!			bc				=> this%boundary%bc_ptr)
	 !
		!time_delay			= 1e-05_dp			
		!time_diff			= 1e-04_dp
		!time_stabilization	= 5e-06_dp			
		!
		!if ( time > (correction+1)*(time_diff) + time_delay) then			
		!			
		!	current_time = time
		!
		!	!# 1D front tracer
		!	do k = cons_inner_loop(3,1),cons_inner_loop(3,2)
		!	do j = cons_inner_loop(2,1),cons_inner_loop(2,2)
		!		max_val = 0.0_dp
		!		do i = cons_inner_loop(1,1),cons_inner_loop(1,2)-1
		!			if(bc%bc_markers(i,j,k) == 0) then	
		!			
		!				!! Grad temp
		!				!if (abs(T%cells(i+1,j,k)-T%cells(i-1,j,k)) > max_val) then
		!				!	max_val = abs(T%cells(i+1,j,k)-T%cells(i-1,j,k))
		!				!	flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
		!				!	flame_front_index = i
		!				!end if
		!			
		!				! max HO2
		!				if (abs(Y%pr(HO2_index)%cells(i,j,k)) > max_val) then
		!					max_val = Y%pr(HO2_index)%cells(i,j,k)
		!					flame_front_coords(j) = (i - 0.5_dp)*cell_size(1) 
		!					flame_front_index = i
		!				end if	
		!			end if
		!		end do
		!	end do
		!	end do
	 !
		!	!left_val	= T%cells(flame_front_index,1,1) - T%cells(flame_front_index-2,1,1)
		!	!right_val	= T%cells(flame_front_index+2,1,1) - T%cells(flame_front_index,1,1)
		!	 
		!	left_val	= Y%pr(HO2_index)%cells(flame_front_index-1,1,1)
		!	right_val	= Y%pr(HO2_index)%cells(flame_front_index+1,1,1)				 
		!	
		!	a = (right_val + left_val - 2.0_dp * max_val)/2.0_dp/cell_size(1)**2
		!	b = (max_val - left_val)/cell_size(1) - a*(2.0_dp*flame_front_coords(1) - cell_size(1))
		!	
		!	current_flame_location = -b/2.0_dp/a
  !
		!	if(correction == 0) then
		!		previous_flame_location = current_flame_location
  !          end if
		!	
  !          boundary = .false.
  !          if(flame_front_index > cons_inner_loop(1,2) - 10) boundary = .true.
  !          
		!	if( (correction /= 0).and.(current_flame_location /=  previous_flame_location) )then 
  !              
		!		flame_velocity = (current_flame_location - previous_flame_location)/(current_time - previous_time)
  !
		!		previous_flame_location = current_flame_location
		!		previous_time = current_time
		!		
  !              previous_av_flame_velocity = sum(flame_velocity_array)
  !              
  !              do i = 19, 1, -1
		!			flame_velocity_array(i+1) = flame_velocity_array(i)  
  !              end do
  !              flame_velocity_array(1) = flame_velocity
  !              
  !              av_flame_velocity = sum(flame_velocity_array)
  !              
		!		if( (correction /= 0).and.(abs(av_flame_velocity - previous_av_flame_velocity) < 1e-03))then 
		!			counter = counter + 1
  !              end if
  !              
  !              write (flame_loc_unit,'(5E14.6)') time, current_flame_location, flame_velocity, av_flame_velocity, abs(av_flame_velocity - previous_av_flame_velocity)
  !
		!	end if
  !
		!	if ((counter > 10).or.boundary) then
		!		stabilized = .true.
		!	end if
		!	
		!	correction = correction + 1
		!	
		!end if	
		!	
		!end associate
  !  end subroutine	
    
    
end module
