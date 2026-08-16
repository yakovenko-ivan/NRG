! CABARET low-Mach projection solver.
!
! CABARET advances the material and momentum equations; a variable-density
! pressure projection enforces the low-Mach divergence constraint on the same
! normal face velocities used by the conservative corrector.  The uniform
! thermodynamic pressure is fixed in open domains and evolves in closed ones.
!
! The material/contact family is reconstructed with temperature and species
! mass fractions, consistently with the validated compressible CABARET solver.
! Face density is then recovered from T_f, Y_f and p_0.  The projection acts
! directly on the face-normal velocities used by the conservative corrector.
! Reduced-gravity buoyancy uses the same ambient-reference convention as FDS.

module cabaret_low_mach_solver_class

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
	use elliptic_multigrid_solver_class

	use benchmarking
	use solver_options_class

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public	:: cabaret_low_mach_solver, cabaret_low_mach_solver_c

	real(dp), parameter :: cabaret_projection_relative_tolerance = 1.0e-4_dp
	real(dp), parameter :: cabaret_projection_absolute_tolerance = 1.0e-10_dp
	real(dp), parameter :: cabaret_projection_density_floor = 1.0e-12_dp
	real(dp), parameter :: cabaret_time_step_growth_limit = 1.25_dp
	real(dp), parameter :: cabaret_closed_pressure_denominator_tolerance = 1.0e-12_dp

	real(dp), parameter :: cabaret_initial_riemann_molar_mass_ratio = 1.01_dp
	real(dp), parameter :: cabaret_initial_riemann_pressure_ratio = 1.05_dp
	real(dp), parameter :: cabaret_initial_riemann_velocity_ratio = 0.05_dp
	real(dp), parameter :: cabaret_initial_contact_velocity_eps = 10.0_dp*tiny(1.0_dp)

	type(field_scalar_cons), target :: p_dyn_field, div_v_field, pressure_correction_field
	type(field_scalar_flow), target :: rho_f_new, p_f_new, e_i_f_new
	type(field_scalar_flow), target :: v_s_f_new, E_f_f_new, T_f_new
	type(field_vector_flow), target :: Y_f_new, v_f_new

#ifdef OMP
	integer(kind=omp_lock_kind), dimension(:,:,:), allocatable :: lock
#endif

	type(timer) :: cabaret_timer
	type(timer) :: cabaret_gas_dynamics_timer
	type(timer) :: cabaret_eos_timer
	type(timer) :: cabaret_chemistry_timer
	type(timer) :: cabaret_diffusion_timer
	type(timer) :: cabaret_heattransfer_timer
	type(timer) :: cabaret_radiation_timer
	type(timer) :: cabaret_viscosity_timer

	type cabaret_low_mach_solver
		logical :: diffusion_flag, viscosity_flag, heat_trans_flag, radiation_flag, reactive_flag
		logical :: CFL_condition_flag, closed_domain

		real(dp) :: courant_fraction, time, time_step
		real(dp) :: thermodynamic_pressure, thermodynamic_pressure_old
		real(dp) :: thermodynamic_pressure_rate_old, thermodynamic_pressure_rate_mid
		real(dp) :: buoyancy_reference_density
		real(dp), dimension(3) :: gravity
		integer :: additional_particles_phases_number

		type(chemical_kinetics_solver) :: chem_kin_solver
		type(diffusion_solver) :: diff_solver
		type(heat_transfer_solver) :: heat_trans_solver
		type(thermal_radiation_solver) :: radiation_solver
		type(viscosity_solver) :: viscosity_solver
		type(table_approximated_real_gas) :: state_eq
		type(dispersed_phase_solver), dimension(:), allocatable :: particles_solver

		type(computational_domain) :: domain
		type(mpi_communications) :: mpi_support
		type(elliptic_multigrid_solver) :: pressure_solver
		type(chemical_properties_pointer) :: chem
		type(thermophysical_properties_pointer) :: thermo
		type(computational_mesh_pointer) :: mesh
		type(boundary_conditions_pointer) :: boundary

		type(field_scalar_cons_pointer) :: rho, T, p, v_s, E_f
		type(field_scalar_cons_pointer) :: mix_mol_mass, mixture_cp, nu, kappa, h_s
		type(field_vector_cons_pointer) :: D, v, Y
		type(field_scalar_cons_pointer) :: p_dyn, div_v, pressure_correction

		type(field_scalar_flow_pointer) :: rho_f_new, p_f_new, e_i_f_new
		type(field_scalar_flow_pointer) :: v_s_f_new, E_f_f_new, T_f_new
		type(field_vector_flow_pointer) :: v_f_new, Y_f_new

		type(field_scalar_cons_pointer) :: E_f_prod_chem, E_f_prod_heat, E_f_prod_rad
		type(field_scalar_cons_pointer) :: E_f_prod_diff
		type(field_vector_cons_pointer) :: v_prod_visc, Y_prod_chem, Y_prod_diff

		type(field_scalar_cons_pointer), dimension(:), allocatable :: rho_prod_particles
		type(field_scalar_cons_pointer), dimension(:), allocatable :: E_f_prod_particles
		type(field_vector_cons_pointer), dimension(:), allocatable :: Y_prod_particles
		type(field_vector_cons_pointer), dimension(:), allocatable :: v_prod_particles

		real(dp), dimension(:,:,:), allocatable :: rho_old
		real(dp), dimension(:,:,:), allocatable :: h_s_old, sensible_enthalpy_advection_old, sensible_enthalpy_source
		real(dp), dimension(:,:,:,:), allocatable :: v_old, Y_old, species_advection_old
		real(dp), dimension(:,:,:,:), allocatable :: momentum_source, species_mass_source, momentum_advection_old
		real(dp), dimension(:,:,:,:), allocatable :: dynamic_pressure_face

		real(dp), dimension(:,:,:,:,:), allocatable :: velocity_face_old, species_face_old
		real(dp), dimension(:,:,:,:), allocatable :: density_face_old

	contains
		procedure, private :: apply_cell_boundary_conditions
		procedure :: solve_problem
		procedure :: calculate_time_step
		procedure :: get_time_step
		procedure :: get_time
		procedure :: set_CFL_coefficient
		procedure, private :: cache_cell_state
		procedure, private :: exchange_conservative_state
		procedure, private :: exchange_face_pressure_density
		procedure, private :: exchange_face_primitive_state
		procedure, private :: exchange_face_thermodynamic_state
		procedure, private :: normalize_face_mass_fractions
		procedure, private :: prepare_gas_dynamics_step
		procedure, private :: calculate_divergence_target
		procedure, private :: thermodynamic_divergence_source
		procedure, private :: reconstruct_projection_provisional_face_velocity
		procedure, private :: apply_low_mach_pressure_projection
		procedure, private :: correct_projection_cell_velocity
		procedure, private :: update_projection_pressure_faces
		procedure, private :: initialize_pressure
		procedure, private :: detect_closed_domain
		procedure, private :: prepare_closed_domain_thermodynamic_pressure_predictor
		procedure, private :: refresh_closed_domain_thermodynamic_pressure_midpoint
		procedure, private :: advance_closed_domain_thermodynamic_pressure_full_step
		procedure, private :: closed_domain_thermodynamic_pressure_rate
		procedure, private :: thermodynamic_pressure_divergence_coefficient
		procedure, private :: evaluate_continuum_source_rates
		procedure, private :: material_characteristic_source_rates
		procedure, private :: predict_material_half_step
		procedure, private :: predict_momentum_half_step
		procedure, private :: correct_material_full_step
		procedure, private :: correct_momentum_full_step
		procedure, private :: temperature_face_value
		procedure, private :: sensible_enthalpy_face_value
		procedure, private :: close_material_face_density
		procedure, private :: refresh_thermodynamic_ghosts
		procedure, private :: update_cell_thermodynamics
		procedure, private :: set_face_thermodynamic_pressure
		procedure, private :: reconstruct_contact_family_face_state
		procedure, private :: finalize_material_face_reconstruction
		procedure, private :: initialize_material_temperature_faces
		procedure, private :: update_flow_thermodynamics
		procedure, private :: finalize_gas_dynamics_step
		procedure, private :: advance_particle_phases
		procedure, private :: cache_face_state_for_next_step
	end type cabaret_low_mach_solver

	interface cabaret_low_mach_solver_c
		module procedure constructor
	end interface

contains
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

	pure real(dp) function limit_quasi_invariant(value, left_value, half_value, right_value, &
		source_shift, lower_clip, upper_clip) result(limited)
		real(dp), intent(in) :: value, left_value, half_value, right_value
		real(dp), intent(in) :: source_shift, lower_clip, upper_clip
		real(dp) :: min_inv, max_inv

		! Source-shifted maximum-principle interval following the
		! balance-characteristic form: the characteristic source term is not
		! interpreted as a spurious new extremum.
		min_inv = min(left_value, half_value, right_value) + source_shift
		max_inv = max(left_value, half_value, right_value) + source_shift
		min_inv = max(min_inv, lower_clip)
		max_inv = min(max_inv, upper_clip)
		if (max_inv < min_inv) then
			max_inv = min_inv
		end if

		limited = min(max(value, min_inv), max_inv)
	end function limit_quasi_invariant

	! Construct field bindings, optional source solvers, persistent CABARET face
	! states, and the low-Mach pressure subsystem.
	type(cabaret_low_mach_solver) function constructor(manager, problem_data_io)
		type(data_manager)						,intent(inout)	:: manager
		type(data_io)							,intent(inout)	:: problem_data_io

		real(dp)								:: calculation_time

		type(field_scalar_cons_pointer)	:: scal_c_ptr
		type(field_vector_cons_pointer)	:: vect_c_ptr
		type(field_tensor_cons_pointer)	:: tens_c_ptr

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

        constructor%gravity                       = manager%solver_options%get_grav_acc()

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
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_molar_mass')
		constructor%mix_mol_mass%s_ptr => scal_c_ptr%s_ptr

		call manager%create_scalar_field(p_dyn_field,'dynamic_pressure','p_dyn')
		constructor%p_dyn%s_ptr => p_dyn_field
		call manager%create_scalar_field(div_v_field, 'velocity_divergence_target', 'div_v')
		constructor%div_v%s_ptr => div_v_field
		call manager%create_scalar_field(pressure_correction_field, &
			'projection_pressure_correction', 'phi')
		constructor%pressure_correction%s_ptr => pressure_correction_field
		constructor%div_v%s_ptr%cells = 0.0_dp
		constructor%pressure_correction%s_ptr%cells = 0.0_dp
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
		constructor%T_f_new%s_ptr 	=> T_f_new

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
			constructor%Y_prod_diff%v_ptr => vect_c_ptr%v_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'diffusivity')
			constructor%D%v_ptr => vect_c_ptr%v_ptr
		end if

		if (constructor%heat_trans_flag .or. &
			constructor%additional_particles_phases_number > 0) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr => scal_c_ptr%s_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'thermal_conductivity')
			constructor%kappa%s_ptr => scal_c_ptr%s_ptr
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
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr => vect_c_ptr%v_ptr
			call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'viscosity')
			constructor%nu%s_ptr => scal_c_ptr%s_ptr
        end if

		constructor%state_eq	=	table_approximated_real_gas_c(manager)
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'sensible_enthalpy')
		constructor%h_s%s_ptr => scal_c_ptr%s_ptr
		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'mixture_cp')
		constructor%mixture_cp%s_ptr => scal_c_ptr%s_ptr

		if(constructor%additional_particles_phases_number /= 0) then
			allocate(constructor%particles_solver(constructor%additional_particles_phases_number))
			allocate(constructor%rho_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%E_f_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%v_prod_particles(constructor%additional_particles_phases_number))
			allocate(constructor%Y_prod_particles(constructor%additional_particles_phases_number))
			do particles_phase_counter = 1, constructor%additional_particles_phases_number
				particles_params = manager%solver_options%get_particles_params(particles_phase_counter)
				constructor%particles_solver(particles_phase_counter) = &
					dispersed_phase_solver_c(manager,particles_params,particles_phase_counter)
				write(var_name,'(A,I2.2)') 'energy_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%E_f_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'density_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%rho_prod_particles(particles_phase_counter)%s_ptr	=> scal_c_ptr%s_ptr
				write(var_name,'(A,I2.2)') 'velocity_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%v_prod_particles(particles_phase_counter)%v_ptr => vect_c_ptr%v_ptr
				write(var_name,'(A,I2.2)') 'concentration_production_particles', particles_phase_counter
				call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,var_name)
				constructor%Y_prod_particles(particles_phase_counter)%v_ptr		=> vect_c_ptr%v_ptr
			end do
		end if

		call manager%get_cons_field_pointer_by_name(scal_c_ptr,vect_c_ptr,tens_c_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_c_ptr%s_ptr

		problem_data_io				= data_io_c(manager,calculation_time)

		if(problem_data_io%get_load_counter() /= 0) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
			call problem_data_io%add_io_scalar_flow_field(constructor%rho_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%p_f_new)
			call problem_data_io%add_io_scalar_flow_field(constructor%E_f_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%Y_f_new)
			call problem_data_io%add_io_vector_flow_field(constructor%v_f_new)
		end if

		call problem_data_io%input_all_data()

		if(problem_data_io%get_load_counter() == 1) then
			call problem_data_io%add_io_scalar_cons_field(constructor%E_f)
			call problem_data_io%add_io_scalar_cons_field(constructor%p_dyn)
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

		allocate(constructor%h_s_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%sensible_enthalpy_advection_old(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))
		allocate(constructor%sensible_enthalpy_source(	cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
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

		allocate(constructor%species_advection_old( species_number, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%momentum_source(	dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%species_mass_source(	species_number					, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%momentum_advection_old(	dimensions						, &
										cons_allocation_bounds(1,1):cons_allocation_bounds(1,2), &
										cons_allocation_bounds(2,1):cons_allocation_bounds(2,2), &
										cons_allocation_bounds(3,1):cons_allocation_bounds(3,2)))

		allocate(constructor%dynamic_pressure_face( dimensions, &
			                                    flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
			                                    flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
			                                    flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		allocate(constructor%density_face_old(		dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		allocate(constructor%velocity_face_old(		dimensions						, &
										dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

		allocate(constructor%species_face_old(		species_number					, &
										dimensions						, &
										flow_allocation_bounds(1,1):flow_allocation_bounds(1,2), &
										flow_allocation_bounds(2,1):flow_allocation_bounds(2,2), &
										flow_allocation_bounds(3,1):flow_allocation_bounds(3,2)))

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

			constructor%density_face_old(:,:,:,:)	= 0.0_dp
			constructor%species_face_old(:,:,:,:,:)	= 0.0_dp
			constructor%velocity_face_old(:,:,:,:,:)	= 0.0_dp
			constructor%momentum_source(:,:,:,:) = 0.0_dp
			constructor%species_mass_source(:,:,:,:) = 0.0_dp
			constructor%momentum_advection_old(:,:,:,:) = 0.0_dp

			do dim = 1, dimensions

				loop = flow_inner_loop

				do dim1 = 1,dimensions
					loop(dim1,2) = flow_inner_loop(dim1,2) - (1 - I_m(dim1,dim))
				end do

				do k = loop(3,1),loop(3,2)
				do j = loop(2,1),loop(2,2)
				do i = loop(1,1),loop(1,2)

							constructor%p_f_new%s_ptr%cells(dim,i,j,k) = 0.5_dp*(constructor%p%s_ptr%cells(i,j,k)+ &
								constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%density_face_old(dim,i,j,k) = 0.5_dp*(constructor%rho%s_ptr%cells(i,j,k)+ &
								constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%E_f_f_new%s_ptr%cells(dim,i,j,k) = 0.5_dp*(constructor%E_f%s_ptr%cells(i,j,k)+ &
								constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							constructor%v_s_f_new%s_ptr%cells(dim,i,j,k) = 0.5_dp*(constructor%v_s%s_ptr%cells(i,j,k)+ &
								constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							spec_summ = 0.0_dp
							do spec = 1, species_number
								constructor%species_face_old(spec,dim,i,j,k) = 0.5_dp*( &
									constructor%Y%v_ptr%pr(spec)%cells(i,j,k)+ &
									constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
								spec_summ = spec_summ + max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)
							end do

							do spec = 1,species_number
								constructor%species_face_old(spec,dim,i,j,k) = &
									max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)/spec_summ
							end do

							do dim1 = 1, dimensions
								constructor%velocity_face_old(dim1,dim,i,j,k) = 0.5_dp*( &
									constructor%v%v_ptr%pr(dim1)%cells(i,j,k)+ &
									constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
							end do

							! initialize the CABARET flux variables by arithmetic averaging.
							! The first predictor step uses these old face states before any
							! characteristic reconstruction is available.  At a H2/air contact,
							! an averaged species_face_old produces a nonphysical mixture face and can
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
								use_riemann_initial_face = &
									molar_mass_ratio > cabaret_initial_riemann_molar_mass_ratio .or. &
									p_ratio > cabaret_initial_riemann_pressure_ratio .or. &
									abs(v_l-v_r) > cabaret_initial_riemann_velocity_ratio*velocity_scale

								if (use_riemann_initial_face) then
									p_floor = 1.0e-12_dp*max(1.0_dp, abs(p_l), abs(p_r))
									rho_floor = 1.0e-12_dp*max(1.0_dp, abs(rho_l), abs(rho_r))
									do spec = 1, species_number
										Y_left_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells( &
										i-I_m(dim,1), j-I_m(dim,2), k-I_m(dim,3))
										Y_right_riemann(spec) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
									end do

									call initial_riemann%set_thermally_perfect_parameters(constructor%thermo%thermo_ptr, &
										rho_l,rho_r,p_l,p_r,v_l,v_r,Y_left_riemann,Y_right_riemann,p_floor,rho_floor)
									call initial_riemann%solve()

									if (initial_riemann%get_success()) then
										u_face_init = initial_riemann%get_velocity()
										u_contact_init = initial_riemann%get_contact_velocity()
										use_left_contact_state = (u_contact_init >= 0.0_dp)
										call initial_riemann%get_mass_fractions(Y_face_riemann)

										constructor%p_f_new%s_ptr%cells(dim,i,j,k)   = initial_riemann%get_pressure()
										constructor%density_face_old(dim,i,j,k) = initial_riemann%get_density()
										constructor%velocity_face_old(dim,dim,i,j,k) = u_face_init
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
											constructor%p_f_new%s_ptr%cells(dim,i,j,k) = p_l
											constructor%density_face_old(dim,i,j,k) = rho_l
											constructor%velocity_face_old(dim,dim,i,j,k) = v_l
										else
											constructor%p_f_new%s_ptr%cells(dim,i,j,k) = p_r
											constructor%density_face_old(dim,i,j,k) = rho_r
											constructor%velocity_face_old(dim,dim,i,j,k) = v_r
										end if
									end if

									if (use_left_contact_state) then
										constructor%E_f_f_new%s_ptr%cells(dim,i,j,k) = &
										constructor%E_f%s_ptr%cells(i-I_m(dim,1), j-I_m(dim,2), k-I_m(dim,3))
										constructor%v_s_f_new%s_ptr%cells(dim,i,j,k) = &
										constructor%v_s%s_ptr%cells(i-I_m(dim,1), j-I_m(dim,2), k-I_m(dim,3))
										do spec = 1, species_number
											if (initial_riemann%get_success()) then
												constructor%species_face_old(spec,dim,i,j,k) = Y_face_riemann(spec)
											else
												constructor%species_face_old(spec,dim,i,j,k) = &
												constructor%Y%v_ptr%pr(spec)%cells( &
												i-I_m(dim,1), j-I_m(dim,2), k-I_m(dim,3))
											end if
										end do
										do dim1 = 1, dimensions
										if (dim1 /= dim) then
											constructor%velocity_face_old(dim1,dim,i,j,k) = &
											constructor%v%v_ptr%pr(dim1)%cells( &
											i-I_m(dim,1), j-I_m(dim,2), k-I_m(dim,3))
										end if
										end do
									else
										constructor%E_f_f_new%s_ptr%cells(dim,i,j,k) = constructor%E_f%s_ptr%cells(i,j,k)
										constructor%v_s_f_new%s_ptr%cells(dim,i,j,k) = constructor%v_s%s_ptr%cells(i,j,k)
										do spec = 1, species_number
											if (initial_riemann%get_success()) then
												constructor%species_face_old(spec,dim,i,j,k) = Y_face_riemann(spec)
											else
												constructor%species_face_old(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
											end if
										end do
										do dim1 = 1, dimensions
										if (dim1 /= dim) then
											constructor%velocity_face_old(dim1,dim,i,j,k) = constructor%v%v_ptr%pr(dim1)%cells(i,j,k)
										end if
										end do
									end if
								end if
							end if

							if (constructor%boundary%bc_ptr%bc_markers(i,j,k) /= 0) then
								constructor%p_f_new%s_ptr%cells(dim,i,j,k)		=	constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								constructor%density_face_old(dim,i,j,k)	=	constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								constructor%E_f_f_new%s_ptr%cells(dim,i,j,k)	=	constructor%E_f%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								constructor%v_s_f_new%s_ptr%cells(dim,i,j,k)	=	constructor%v_s%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))

								spec_summ = 0.0_dp
								do spec = 1, species_number
									constructor%species_face_old(spec,dim,i,j,k) =  constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
									spec_summ = spec_summ + max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)
								end do

								do spec = 1,species_number
									constructor%species_face_old(spec,dim,i,j,k) = &
									max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)/spec_summ
								end do

								do dim1 = 1, dimensions
									constructor%velocity_face_old(dim1,dim,i,j,k) =	constructor%v%v_ptr%pr(dim1)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
								end do
							end if

							if (constructor%boundary%bc_ptr%bc_markers(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) /= 0) then
								constructor%p_f_new%s_ptr%cells(dim,i,j,k)		=	constructor%p%s_ptr%cells(i,j,k)
								constructor%density_face_old(dim,i,j,k)	=	constructor%rho%s_ptr%cells(i,j,k)
								constructor%E_f_f_new%s_ptr%cells(dim,i,j,k)	=	constructor%E_f%s_ptr%cells(i,j,k)
								constructor%v_s_f_new%s_ptr%cells(dim,i,j,k)	=	constructor%v_s%s_ptr%cells(i,j,k)

								spec_summ = 0.0_dp
								do spec = 1, species_number
									constructor%species_face_old(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
									spec_summ = spec_summ + max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)
								end do

								do spec = 1,species_number
									constructor%species_face_old(spec,dim,i,j,k) = &
									max(constructor%species_face_old(spec,dim,i,j,k), 0.0_dp)/spec_summ
								end do

								do dim1 = 1, dimensions
									constructor%velocity_face_old(dim1,dim,i,j,k) =	constructor%v%v_ptr%pr(dim1)%cells(i,j,k)
								end do

							end if
				end do
				end do
				end do
			end do

			do spec = 1,species_number
				constructor%Y_f_new%v_ptr%pr(spec)%cells(:,:,:,:) = constructor%species_face_old(spec,:,:,:,:)
			end do

			do dim = 1,dimensions
				constructor%v_f_new%v_ptr%pr(dim)%cells(:,:,:,:) = constructor%velocity_face_old(dim,:,:,:,:)
			end do

			constructor%rho_f_new%s_ptr%cells(:,:,:,:)	= constructor%density_face_old

		end if

		call constructor%mpi_support%exchange_flow_scalar_field(constructor%p_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%rho_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%E_f_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_scalar_field(constructor%v_s_f_new%s_ptr)
		call constructor%mpi_support%exchange_flow_vector_field(constructor%Y_f_new%v_ptr)
		call constructor%mpi_support%exchange_flow_vector_field(constructor%v_f_new%v_ptr)

		do spec = 1,species_number
			constructor%species_face_old(spec,:,:,:,:)	= constructor%Y_f_new%v_ptr%pr(spec)%cells
		end do

		do dim = 1,dimensions
			constructor%velocity_face_old(dim,:,:,:,:)     = constructor%v_f_new%v_ptr%pr(dim)%cells
		end do

		constructor%time		        =   calculation_time
		constructor%time_step	        =   manager%solver_options%get_initial_time_step()

		constructor%thermodynamic_pressure = 0.0_dp
		constructor%thermodynamic_pressure_old = 0.0_dp
		constructor%thermodynamic_pressure_rate_old = 0.0_dp
		constructor%thermodynamic_pressure_rate_mid = 0.0_dp
		constructor%closed_domain = .false.
		constructor%dynamic_pressure_face = 0.0_dp
		constructor%sensible_enthalpy_advection_old = 0.0_dp
		constructor%sensible_enthalpy_source = 0.0_dp
		constructor%species_advection_old = 0.0_dp

		call constructor%state_eq%apply_state_equation_flow_variables_for_IC()
		call constructor%apply_cell_boundary_conditions()

		! FDS-compatible reference density for the reduced-gravity force.  The
		! upper outer fluid cell is normally the undisturbed ambient state.
		constructor%buoyancy_reference_density = constructor%rho%s_ptr%cells(cons_inner_loop(1,2), &
			cons_inner_loop(2,2),cons_inner_loop(3,2))
		if (constructor%buoyancy_reference_density <= cabaret_projection_density_floor) &
			error stop 'CABARET projection: non-positive buoyancy reference density'

		call constructor%initialize_pressure(problem_data_io%get_load_counter() == 1)
        constructor%density_face_old   = constructor%rho_f_new%s_ptr%cells

        call manager%create_timer(cabaret_timer                 ,'CABARET solver time'              , 'sol_t')
        call manager%create_timer(cabaret_gas_dynamics_timer    ,'CABARET gas dynamics time'        , 'gd_t')
        call manager%create_timer(cabaret_eos_timer             ,'CABARET eos solver time'          , 'eos_t')
        call manager%create_timer(cabaret_chemistry_timer       ,'CABARET chemistry solver time'    , 'chem_t')
        call manager%create_timer(cabaret_diffusion_timer       ,'CABARET diffusion solver time'    , 'diff_t')
        call manager%create_timer(cabaret_heattransfer_timer    ,'CABARET heattransfer solver time' , 'ht_t')
        call manager%create_timer(cabaret_radiation_timer       ,'CABARET radiation solver time'    , 'rad_t')
        call manager%create_timer(cabaret_viscosity_timer       ,'CABARET viscosity solver time'    , 'visc_t')
	end function

! Advance one complete time step: explicit CABARET prediction, two face
! projections around material reconstruction, conservative correction, and
! optional particle coupling.
subroutine solve_problem(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    call cabaret_timer%tic()

    call this%prepare_gas_dynamics_step()
    call this%evaluate_continuum_source_rates()
    call this%prepare_closed_domain_thermodynamic_pressure_predictor()

    call cabaret_gas_dynamics_timer%tic()
    call this%predict_material_half_step()
    call this%update_cell_thermodynamics()
    call this%predict_momentum_half_step()
    call this%update_cell_thermodynamics()
    call this%refresh_closed_domain_thermodynamic_pressure_midpoint()
    call this%exchange_conservative_state()
    call this%apply_cell_boundary_conditions()
    call cabaret_gas_dynamics_timer%toc()

    call cabaret_gas_dynamics_timer%tic()
    call this%reconstruct_projection_provisional_face_velocity()
    call this%calculate_divergence_target()
    call this%apply_low_mach_pressure_projection(.false.)
    call this%reconstruct_contact_family_face_state()
    call this%finalize_material_face_reconstruction()
    call this%close_material_face_density()
    call this%update_flow_thermodynamics()
    call this%apply_low_mach_pressure_projection(.true.)
    call cabaret_gas_dynamics_timer%toc()

    call cabaret_gas_dynamics_timer%tic()
    call this%advance_closed_domain_thermodynamic_pressure_full_step()
    call this%correct_material_full_step()
    call this%update_cell_thermodynamics()
    call this%correct_momentum_full_step()
    call this%correct_projection_cell_velocity()
    call this%update_cell_thermodynamics()
    call this%finalize_gas_dynamics_step()
    call cabaret_gas_dynamics_timer%toc()

    ! Dispersed-phase rates were included in the predictor/corrector source
    ! assembly, so no post-step direct state increment is required.

    call cabaret_gas_dynamics_timer%tic()
    call this%cache_face_state_for_next_step()
    if (this%CFL_condition_flag) call this%calculate_time_step()
    this%time = this%time + this%time_step
    call cabaret_gas_dynamics_timer%toc(new_iter=.true.)
    call cabaret_timer%toc(new_iter=.true.)
end subroutine solve_problem

	! CABARET gas-dynamics stages

subroutine prepare_gas_dynamics_step(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    call this%apply_cell_boundary_conditions()
    call this%refresh_thermodynamic_ghosts()
    call this%exchange_conservative_state()
    call this%cache_cell_state()
end subroutine prepare_gas_dynamics_step

! Split pressure into a spatially uniform thermodynamic component and a
! dynamic projection component.  Restarted runs preserve the latter.
subroutine initialize_pressure(this,reset_dynamic_pressure)
#ifdef mpi
    use MPI
#endif
    class(cabaret_low_mach_solver), intent(inout) :: this
    logical, intent(in) :: reset_dynamic_pressure

    integer :: i, j, k, dim
    integer, dimension(3,2) :: cell_loop, face_loop
    real(dp) :: pressure_sum_local, pressure_sum_global
    real(dp) :: cell_count_local, cell_count_global
#ifdef mpi
    integer :: error, mpi_communicator
#endif

    cell_loop = this%domain%get_local_inner_cells_bounds()
    face_loop = this%domain%get_local_utter_faces_bounds()
    pressure_sum_local = 0.0_dp
    cell_count_local = 0.0_dp

    associate(p => this%p%s_ptr, bc => this%boundary%bc_ptr)
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (bc%bc_markers(i,j,k) == 0) then
                        pressure_sum_local = pressure_sum_local+p%cells(i,j,k)
                        cell_count_local = cell_count_local + 1.0_dp
                    end if
                end do
            end do
        end do
    end associate

#ifdef mpi
    mpi_communicator = this%domain%get_mpi_communicator()
    call mpi_allreduce(pressure_sum_local,pressure_sum_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_communicator,error)
    call mpi_allreduce(cell_count_local,cell_count_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_communicator,error)
#else
    pressure_sum_global = pressure_sum_local
    cell_count_global = cell_count_local
#endif
    if (cell_count_global <= 0.0_dp) error stop 'CABARET projection: empty domain'

    this%thermodynamic_pressure = pressure_sum_global/cell_count_global
    if (this%thermodynamic_pressure <= 0.0_dp) error stop 'CABARET projection: non-positive thermodynamic pressure'
    this%thermodynamic_pressure_old = this%thermodynamic_pressure
    this%thermodynamic_pressure_rate_old = 0.0_dp
    this%thermodynamic_pressure_rate_mid = 0.0_dp
    this%dynamic_pressure_face = 0.0_dp
    if (reset_dynamic_pressure) then
        this%p_dyn%s_ptr%cells = 0.0_dp
        this%p_f_new%s_ptr%cells = this%thermodynamic_pressure
    else
        call this%mpi_support%exchange_conservative_scalar_field(this%p_dyn%s_ptr)
        call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
    end if

    do dim = 1, this%domain%get_domain_dimensions()
        do k = face_loop(3,1), face_loop(3,2)
            do j = face_loop(2,1), face_loop(2,2)
                do i = face_loop(1,1), face_loop(1,2)
                    this%dynamic_pressure_face(dim,i,j,k) = this%p_f_new%s_ptr%cells(dim,i,j,k)-this%thermodynamic_pressure
                end do
            end do
        end do
    end do
    call this%detect_closed_domain()
end subroutine initialize_pressure

! A domain is thermodynamically closed only when every physical boundary is
! a wall or symmetry plane; open boundaries keep p_0 fixed.
subroutine detect_closed_domain(this)
#ifdef mpi
    use MPI
#endif
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, i, j, k, dim, side, sign, boundary_number
    integer :: gi, gj, gk, closed_local, closed_global
    integer, dimension(3,2) :: loop
    character(len=20) :: boundary_name
#ifdef mpi
    integer :: error, mpi_communicator
#endif

    dimensions = this%domain%get_domain_dimensions()
    loop = this%domain%get_local_inner_cells_bounds()
    closed_local = 1
    associate(bc => this%boundary%bc_ptr)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do dim = 1, dimensions
                        do side = 1, 2
                            sign = (-1)**side
                            gi = i+sign*I_m(dim,1)
                            gj = j+sign*I_m(dim,2)
                            gk = k+sign*I_m(dim,3)
                            boundary_number = bc%bc_markers(gi,gj,gk)
                            if (boundary_number == 0) cycle
                            boundary_name = bc%boundary_types(boundary_number)%get_type_name()
                            if (boundary_name /= 'wall' .and. boundary_name /= 'symmetry_plane') closed_local = 0
                        end do
                    end do
                end do
            end do
        end do
    end associate
#ifdef mpi
    mpi_communicator = this%domain%get_mpi_communicator()
    call mpi_allreduce(closed_local,closed_global,1,MPI_INTEGER,MPI_MIN,mpi_communicator,error)
#else
    closed_global = closed_local
#endif
    this%closed_domain = (closed_global == 1)
end subroutine detect_closed_domain

	real(dp) function thermodynamic_pressure_divergence_coefficient(this, rho_state, p_state, Y_state) result(beta)
		class(cabaret_low_mach_solver), intent(in) :: this
		real(dp), intent(in) :: rho_state, p_state
		real(dp), dimension(:), intent(in) :: Y_state
		real(dp) :: mixture_molar_mass, temperature_state, cp_molar, cp_mass
		real(dp), dimension(size(Y_state)) :: X_mole

		if (rho_state <= 0.0_dp .or. p_state <= 0.0_dp) &
			error stop 'CABARET projection: invalid thermodynamic-pressure state'
		mixture_molar_mass = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_state)
		temperature_state = p_state*mixture_molar_mass/(rho_state*r_gase_J)
		call this%thermo%thermo_ptr%mole_fractions_from_mass_fractions(Y_state, X_mole)
		cp_molar = this%thermo%thermo_ptr%mixture_cp_molar(temperature_state, X_mole)
		cp_mass = cp_molar/mixture_molar_mass
		if (temperature_state <= 0.0_dp .or. cp_mass <= 0.0_dp) &
			error stop 'CABARET projection: invalid heat capacity in pressure constraint'
		beta = 1.0_dp/p_state - 1.0_dp/(rho_state*cp_mass*temperature_state)
	end function thermodynamic_pressure_divergence_coefficient

! Enforce the volume-integrated divergence compatibility condition in a
! closed domain to obtain dp_0/dt.
real(dp) function closed_domain_thermodynamic_pressure_rate(this,old_state) result(rate)
#ifdef mpi
    use MPI
#endif
    class(cabaret_low_mach_solver), intent(in) :: this
    logical, intent(in) :: old_state

    integer :: species_number, i, j, k, spec, geometry_power
    integer, dimension(3,2) :: loop
    real(dp) :: rho_state, temperature_state, target, beta, volume, base_volume, radius
    real(dp) :: source_integral_local, coefficient_integral_local, volume_local
    real(dp) :: source_integral, coefficient_integral, total_volume
    real(dp) :: denominator_scale
    real(dp), dimension(:), allocatable :: Y_state
#ifdef mpi
    integer :: error, mpi_communicator
#endif

    rate = 0.0_dp
    if (.not. this%closed_domain) return

    species_number = this%chem%chem_ptr%species_number
    loop = this%domain%get_local_inner_cells_bounds()
    base_volume = this%mesh%mesh_ptr%get_cell_volume()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        geometry_power = 0
    case('cylindrical')
        geometry_power = 1
    case('spherical')
        geometry_power = 2
    case default
        geometry_power = 0
    end select

    allocate(Y_state(species_number))
    source_integral_local = 0.0_dp
    coefficient_integral_local = 0.0_dp
    volume_local = 0.0_dp
    associate(rho => this%rho%s_ptr, T => this%T%s_ptr, Y => this%Y%v_ptr, bc => this%boundary%bc_ptr)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (old_state) then
                        rho_state = this%rho_old(i,j,k)
                        do spec = 1, species_number
                            Y_state(spec) = this%Y_old(spec,i,j,k)
                        end do
                    else
                        rho_state = rho%cells(i,j,k)
                        do spec = 1, species_number
                            Y_state(spec) = Y%pr(spec)%cells(i,j,k)
                        end do
                    end if
                    temperature_state = T%cells(i,j,k)
                    call thermodynamic_divergence_source(this,rho_state,temperature_state,Y_state, &
                        this%sensible_enthalpy_source(i,j,k),this%species_mass_source(:,i,j,k),target)
                    beta = this%thermodynamic_pressure_divergence_coefficient(rho_state,this%thermodynamic_pressure,Y_state)
                    volume = base_volume
                    if (geometry_power > 0) then
                        radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
                        volume = volume*radius**geometry_power
                    end if
                    source_integral_local = source_integral_local+target*volume
                    coefficient_integral_local = coefficient_integral_local+beta*volume
                    volume_local = volume_local+volume
                end do
            end do
        end do
    end associate
    deallocate(Y_state)
#ifdef mpi
    mpi_communicator = this%domain%get_mpi_communicator()
    call mpi_allreduce(source_integral_local,source_integral,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_communicator,error)
    call mpi_allreduce(coefficient_integral_local,coefficient_integral,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_communicator,error)
    call mpi_allreduce(volume_local,total_volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_communicator,error)
#else
    source_integral = source_integral_local
    coefficient_integral = coefficient_integral_local
    total_volume = volume_local
#endif
    denominator_scale = max(total_volume/max(this%thermodynamic_pressure,tiny(1.0_dp)),tiny(1.0_dp))
    if (abs(coefficient_integral) <= cabaret_closed_pressure_denominator_tolerance*denominator_scale) &
        error stop 'CABARET projection: singular closed-domain pressure constraint'
    rate = source_integral/coefficient_integral
end function closed_domain_thermodynamic_pressure_rate

	subroutine prepare_closed_domain_thermodynamic_pressure_predictor(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		this%thermodynamic_pressure_old = this%thermodynamic_pressure
		this%thermodynamic_pressure_rate_old = 0.0_dp
		this%thermodynamic_pressure_rate_mid = 0.0_dp
		if (.not. this%closed_domain) return

		this%thermodynamic_pressure_rate_old = &
			this%closed_domain_thermodynamic_pressure_rate(.true.)
		this%thermodynamic_pressure_rate_mid = &
			this%thermodynamic_pressure_rate_old
		this%thermodynamic_pressure = this%thermodynamic_pressure_old + &
			0.5_dp*this%time_step*this%thermodynamic_pressure_rate_old
		if (this%thermodynamic_pressure <= 0.0_dp) &
			error stop 'CABARET projection: non-positive predicted thermodynamic pressure'
	end subroutine prepare_closed_domain_thermodynamic_pressure_predictor

	subroutine refresh_closed_domain_thermodynamic_pressure_midpoint(this)
		class(cabaret_low_mach_solver), intent(inout) :: this
		if (.not. this%closed_domain) return
		this%thermodynamic_pressure_rate_mid = &
			this%closed_domain_thermodynamic_pressure_rate(.false.)
	end subroutine refresh_closed_domain_thermodynamic_pressure_midpoint

	subroutine advance_closed_domain_thermodynamic_pressure_full_step(this)
		class(cabaret_low_mach_solver), intent(inout) :: this
		if (.not. this%closed_domain) return
		this%thermodynamic_pressure = this%thermodynamic_pressure_old + &
			this%time_step*this%thermodynamic_pressure_rate_mid
		if (this%thermodynamic_pressure <= 0.0_dp) &
			error stop 'CABARET projection: non-positive corrected thermodynamic pressure'
	end subroutine advance_closed_domain_thermodynamic_pressure_full_step

subroutine thermodynamic_divergence_source(this,rho_state,temperature,Y_state,rhoHs_source,rhoY_source,target)
    class(cabaret_low_mach_solver), intent(in) :: this
    real(dp), intent(in) :: rho_state, temperature, rhoHs_source
    real(dp), dimension(:), intent(in) :: Y_state, rhoY_source
    real(dp), intent(out) :: target
    integer :: species_number, spec
    real(dp) :: molar_mass, cp_molar, cp_mass, h_species
    real(dp), dimension(size(Y_state)) :: X_state

    species_number = size(Y_state)
    molar_mass = this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_state)
    call this%thermo%thermo_ptr%mole_fractions_from_mass_fractions(Y_state,X_state)
    cp_molar = this%thermo%thermo_ptr%mixture_cp_molar(temperature,X_state)
    cp_mass = cp_molar/molar_mass
    if (rho_state <= 0.0_dp .or. temperature <= 0.0_dp .or. cp_mass <= 0.0_dp) &
        error stop 'CABARET projection: invalid state in thermodynamic divergence source'
    target = rhoHs_source/(rho_state*cp_mass*temperature)
    do spec = 1, species_number
        h_species = (this%thermo%thermo_ptr%specie_enthalpy_molar(temperature,spec)- &
                     this%thermo%thermo_ptr%specie_enthalpy_molar(T_ref,spec))/ &
                     this%thermo%thermo_ptr%molar_masses(spec)
        target = target+(molar_mass/(rho_state*this%thermo%thermo_ptr%molar_masses(spec))- &
                 h_species/(rho_state*cp_mass*temperature))*rhoY_source(spec)
    end do
end subroutine thermodynamic_divergence_source

! Build the low-Mach dilatation target from sensible-enthalpy and species
! sources, including the closed-domain thermodynamic-pressure correction.
subroutine calculate_divergence_target(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: species_number, i, j, k, spec
    integer, dimension(3,2) :: loop
    real(dp) :: target_value
    real(dp), dimension(:), allocatable :: Y_cell

    species_number = this%chem%chem_ptr%species_number
    loop = this%domain%get_local_inner_cells_bounds()
    associate(rho => this%rho%s_ptr,T => this%T%s_ptr,Y => this%Y%v_ptr, &
            target => this%div_v%s_ptr,bc => this%boundary%bc_ptr)
        !$omp parallel default(shared) private(i,j,k,spec,target_value,Y_cell)
        allocate(Y_cell(species_number))
        !$omp do collapse(3) schedule(static)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    do spec = 1, species_number
                        Y_cell(spec) = Y%pr(spec)%cells(i,j,k)
                    end do
                    call this%thermodynamic_divergence_source(rho%cells(i,j,k),T%cells(i,j,k),Y_cell, &
                        this%sensible_enthalpy_source(i,j,k),this%species_mass_source(:,i,j,k),target_value)
                    if (this%closed_domain) target_value = target_value- &
                        this%thermodynamic_pressure_divergence_coefficient(rho%cells(i,j,k),this%thermodynamic_pressure,Y_cell)* &
                        this%thermodynamic_pressure_rate_mid
                    target%cells(i,j,k) = target_value
                end do
            end do
        end do
        !$omp end do
        deallocate(Y_cell)
        !$omp end parallel
    end associate
    call this%mpi_support%exchange_conservative_scalar_field(this%div_v%s_ptr)
end subroutine calculate_divergence_target

	! CABARET pressure-projection subsystem
! Extrapolate CABARET face velocities from half-step cell values, apply the
! source-shifted maximum principle, and impose physical face velocities.
subroutine reconstruct_projection_provisional_face_velocity(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: dimensions, species_number, i,j,k,dim,comp,spec,plus,sign,bound_number
    integer :: fi,fj,fk,gi,gj,gk
    integer, dimension(3,2) :: cons_inner_loop, cons_utter_loop, loop
    real(dp) :: u_lower,u_upper,u_normal_lower,u_normal_upper,u_half,u_old_cell
    real(dp) :: candidate_lower,candidate_upper,lower_bound,upper_bound,source_shift
    real(dp) :: farfield_velocity,species_sum
    character(len=20) :: boundary_type_name

    dimensions=this%domain%get_domain_dimensions()
    species_number=this%chem%chem_ptr%species_number
    cons_inner_loop=this%domain%get_local_inner_cells_bounds()
    cons_utter_loop=this%domain%get_local_utter_cells_bounds()

    ! Preserve the CABARET old-time flux variables as the default state.  The
    ! one-cell extrapolation below changes only faces selected by the sign of
    ! the old normal face velocity; only the CABARET flux variables enter.
    this%rho_f_new%s_ptr%cells=this%density_face_old
    this%p_f_new%s_ptr%cells=this%thermodynamic_pressure
    do comp=1,dimensions
        this%v_f_new%v_ptr%pr(comp)%cells=this%velocity_face_old(comp,:,:,:,:)
    end do
    do spec=1,species_number
        this%Y_f_new%v_ptr%pr(spec)%cells=this%species_face_old(spec,:,:,:,:)
    end do

    associate(v=>this%v%v_ptr,Y=>this%Y%v_ptr,rho=>this%rho%s_ptr,bc=>this%boundary%bc_ptr)
        do dim=1,dimensions
            loop(3,1)=cons_utter_loop(3,1)*I_m(dim,3)+cons_inner_loop(3,1)*(1-I_m(dim,3))
            loop(3,2)=cons_utter_loop(3,2)*I_m(dim,3)+cons_inner_loop(3,2)*(1-I_m(dim,3))
            loop(2,1)=cons_utter_loop(2,1)*I_m(dim,2)+cons_inner_loop(2,1)*(1-I_m(dim,2))
            loop(2,2)=cons_utter_loop(2,2)*I_m(dim,2)+cons_inner_loop(2,2)*(1-I_m(dim,2))
            loop(1,1)=cons_utter_loop(1,1)*I_m(dim,1)+cons_inner_loop(1,1)*(1-I_m(dim,1))
            loop(1,2)=cons_utter_loop(1,2)*I_m(dim,1)+cons_inner_loop(1,2)*(1-I_m(dim,1))
            do k=loop(3,1),loop(3,2)
                do j=loop(2,1),loop(2,2)
                    do i=loop(1,1),loop(1,2)
                        if (bc%bc_markers(i,j,k)/=0) cycle
                        u_normal_lower=this%velocity_face_old(dim,dim,i,j,k)
                        u_normal_upper=this%velocity_face_old(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        do comp=1,dimensions
                            u_lower=this%velocity_face_old(comp,dim,i,j,k)
                            u_upper=this%velocity_face_old(comp,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                            u_half=v%pr(comp)%cells(i,j,k)
                            u_old_cell=this%v_old(comp,i,j,k)
                            candidate_lower=2.0_dp*u_half-u_upper
                            candidate_upper=2.0_dp*u_half-u_lower
                            ! The half-step cell velocity already contains one half of every
                            ! momentum source.  Shift the CABARET admissible interval by the
                            ! corresponding full-step acceleration so uniform buoyant acceleration
                            ! is not clipped back to the old face velocity.
                            source_shift=this%time_step*this%momentum_source(comp,i,j,k)/ &
                                max(this%rho_old(i,j,k),cabaret_projection_density_floor)
                            lower_bound=min(u_lower,u_old_cell,u_upper)+source_shift
                            upper_bound=max(u_lower,u_old_cell,u_upper)+source_shift
                            candidate_lower=min(upper_bound,max(lower_bound,candidate_lower))
                            candidate_upper=min(upper_bound,max(lower_bound,candidate_upper))
                            if ((I_m(dim,1)*i+I_m(dim,2)*j+I_m(dim,3)*k)/=cons_utter_loop(dim,1)) then
                                gi=i-I_m(dim,1)
                                gj=j-I_m(dim,2)
                                gk=k-I_m(dim,3)
                                if (bc%bc_markers(gi,gj,gk)==0 .and. u_normal_lower<0.0_dp) &
                                    this%v_f_new%v_ptr%pr(comp)%cells(dim,i,j,k)=candidate_lower
                            end if
                            if ((I_m(dim,1)*i+I_m(dim,2)*j+I_m(dim,3)*k)/=cons_utter_loop(dim,2)) then
                                gi=i+I_m(dim,1)
                                gj=j+I_m(dim,2)
                                gk=k+I_m(dim,3)
                                if (bc%bc_markers(gi,gj,gk)==0 .and. u_normal_upper>=0.0_dp) &
                                    this%v_f_new%v_ptr%pr(comp)%cells(dim,gi,gj,gk)=candidate_upper
                            end if
                        end do
                    end do
                end do
            end do
        end do

        ! Physical boundary faces.  Solid/symmetry and prescribed-velocity inlet
        ! faces are Neumann pressure boundaries; a free outlet is the pressure
        ! Dirichlet boundary in the projection solve.
        do k=cons_inner_loop(3,1),cons_inner_loop(3,2)
            do j=cons_inner_loop(2,1),cons_inner_loop(2,2)
                do i=cons_inner_loop(1,1),cons_inner_loop(1,2)
                    if (bc%bc_markers(i,j,k)/=0) cycle
                    do dim=1,dimensions
                        do plus=1,2
                            sign=(-1)**plus
                            gi=i+sign*I_m(dim,1)
                            gj=j+sign*I_m(dim,2)
                            gk=k+sign*I_m(dim,3)
                            bound_number=bc%bc_markers(gi,gj,gk)
                            if (bound_number==0) cycle
                            if (sign<0) then
                                fi=i
                                fj=j
                                fk=k
                            else
                                fi=i+I_m(dim,1)
                                fj=j+I_m(dim,2)
                                fk=k+I_m(dim,3)
                            end if
                            boundary_type_name=bc%boundary_types(bound_number)%get_type_name()
                            do comp=1,dimensions
                                this%v_f_new%v_ptr%pr(comp)%cells(dim,fi,fj,fk)=v%pr(comp)%cells(i,j,k)
                            end do
                            select case(boundary_type_name)
                            case('wall','symmetry_plane')
                                this%v_f_new%v_ptr%pr(dim)%cells(dim,fi,fj,fk)=0.0_dp
                            case('inlet')
                                farfield_velocity=bc%boundary_types(bound_number)%get_farfield_velocity()
                                this%v_f_new%v_ptr%pr(dim)%cells(dim,fi,fj,fk)=farfield_velocity
                            case('outlet')
                                ! zero-gradient provisional velocity; pressure is imposed by phi=0
                            case default
                                error stop 'CABARET projection: unsupported boundary type'
                            end select
                            this%rho_f_new%s_ptr%cells(dim,fi,fj,fk)=max(rho%cells(i,j,k),cabaret_projection_density_floor)
                            this%p_f_new%s_ptr%cells(dim,fi,fj,fk)=this%thermodynamic_pressure
                            species_sum=0.0_dp
                            do spec=1,species_number
                                this%Y_f_new%v_ptr%pr(spec)%cells(dim,fi,fj,fk)=max(Y%pr(spec)%cells(i,j,k),0.0_dp)
                                species_sum=species_sum+this%Y_f_new%v_ptr%pr(spec)%cells(dim,fi,fj,fk)
                            end do
                            if (species_sum>0.0_dp) then
                                do spec=1,species_number
                                    this%Y_f_new%v_ptr%pr(spec)%cells(dim,fi,fj,fk)= &
                                        this%Y_f_new%v_ptr%pr(spec)%cells(dim,fi,fj,fk)/species_sum
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end associate
    call this%exchange_face_primitive_state()
end subroutine reconstruct_projection_provisional_face_velocity

subroutine apply_low_mach_pressure_projection(this,accumulate_pressure)
    class(cabaret_low_mach_solver), intent(inout) :: this
    logical, intent(in) :: accumulate_pressure
    type(elliptic_solver_statistics) :: statistics
    type(elliptic_boundary_data) :: elliptic_boundary
    integer :: dimensions,geometry_power,i,j,k,dim,dim1,ii,jj,kk,il,jl,kl,ir,jr,kr
    integer :: bound_l,bound_r,boundary_marker
    integer, dimension(3,2) :: loop,flow_loop,face_loop
    real(dp), dimension(3) :: dx
    real(dp), allocatable :: solution(:,:,:),rhs(:,:,:),cell_volume(:,:,:)
    real(dp), allocatable :: conductance_x(:,:,:),conductance_y(:,:,:),conductance_z(:,:,:)
    logical, allocatable :: active(:,:,:)
    real(dp) :: base_volume,div_u,rho_face,grad_phi,beta_face
    real(dp) :: radius,radius_low,radius_high,radius_face,metric_center,metric_low,metric_high,metric_face
    real(dp) :: lower_face_velocity,upper_face_velocity
    character(len=20) :: type_l,type_r,boundary_name

    dimensions=this%domain%get_domain_dimensions()
    loop=this%domain%get_local_inner_cells_bounds()
    flow_loop=this%domain%get_local_inner_faces_bounds()
    dx=this%mesh%mesh_ptr%get_cell_edges_length()
    base_volume=this%mesh%mesh_ptr%get_cell_volume()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        geometry_power=0
    case('cylindrical')
        geometry_power=1
    case('spherical')
        geometry_power=2
    case default
        error stop 'CABARET projection: unsupported coordinate system'
    end select

    allocate(solution(loop(1,2)-loop(1,1)+1,loop(2,2)-loop(2,1)+1,loop(3,2)-loop(3,1)+1))
    allocate(rhs,source=solution)
    allocate(cell_volume,source=solution)
    allocate(active(size(solution,1),size(solution,2),size(solution,3)))
    allocate(conductance_x(size(solution,1)+1,size(solution,2),size(solution,3)))
    allocate(conductance_y(size(solution,1),size(solution,2)+1,size(solution,3)))
    allocate(conductance_z(size(solution,1),size(solution,2),size(solution,3)+1))
    allocate(elliptic_boundary%type_x(size(conductance_x,1),size(conductance_x,2),size(conductance_x,3)))
    allocate(elliptic_boundary%value_x(size(conductance_x,1),size(conductance_x,2),size(conductance_x,3)))
    allocate(elliptic_boundary%type_y(size(conductance_y,1),size(conductance_y,2),size(conductance_y,3)))
    allocate(elliptic_boundary%value_y(size(conductance_y,1),size(conductance_y,2),size(conductance_y,3)))
    allocate(elliptic_boundary%type_z(size(conductance_z,1),size(conductance_z,2),size(conductance_z,3)))
    allocate(elliptic_boundary%value_z(size(conductance_z,1),size(conductance_z,2),size(conductance_z,3)))
    solution=0.0_dp
    rhs=0.0_dp
    cell_volume=base_volume
    active=.false.
    conductance_x=0.0_dp
    conductance_y=0.0_dp
    conductance_z=0.0_dp
    elliptic_boundary%type_x=elliptic_bc_internal
    elliptic_boundary%type_y=elliptic_bc_internal
    elliptic_boundary%type_z=elliptic_bc_internal
    elliptic_boundary%value_x=0.0_dp
    elliptic_boundary%value_y=0.0_dp
    elliptic_boundary%value_z=0.0_dp

    call this%exchange_face_pressure_density()
    call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
    associate(bc=>this%boundary%bc_ptr,target=>this%div_v%s_ptr,density_face_old=>this%rho_f_new%s_ptr, &
            mesh=>this%mesh%mesh_ptr)
        do k=loop(3,1),loop(3,2)
            do j=loop(2,1),loop(2,2)
                do i=loop(1,1),loop(1,2)
                    ii=i-loop(1,1)+1
                    jj=j-loop(2,1)+1
                    kk=k-loop(3,1)+1
                    if(bc%bc_markers(i,j,k)/=0)cycle
                    active(ii,jj,kk)=.true.
                    radius=mesh%mesh(1,i,j,k)
                    metric_center=1.0_dp
                    if(geometry_power>0)metric_center=max(radius,0.0_dp)**geometry_power
                    cell_volume(ii,jj,kk)=base_volume*metric_center
                    div_u=0.0_dp
                    do dim=1,dimensions
                        lower_face_velocity=this%v_f_new%v_ptr%pr(dim)%cells(dim,i,j,k)
                        upper_face_velocity=this%v_f_new%v_ptr%pr(dim)%cells(dim, &
                            i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        if(dim==1.and.geometry_power>0)then
                            radius_low=max(radius-0.5_dp*dx(1),0.0_dp)
                            radius_high=radius+0.5_dp*dx(1)
                            metric_low=radius_low**geometry_power
                            metric_high=radius_high**geometry_power
                            div_u=div_u+(metric_high*upper_face_velocity-metric_low*lower_face_velocity)/ &
                                max(metric_center*dx(1),tiny(1.0_dp))
                        else
                            div_u=div_u+(upper_face_velocity-lower_face_velocity)/dx(dim)
                        end if
                    end do
                    rhs(ii,jj,kk)=target%cells(i,j,k)-div_u
                end do
            end do
        end do

        do kk=1,size(solution,3)
            do jj=1,size(solution,2)
                do ii=1,size(solution,1)+1
                    i=loop(1,1)+ii-1
                    j=loop(2,1)+jj-1
                    k=loop(3,1)+kk-1
                    rho_face=max(density_face_old%cells(1,i,j,k),cabaret_projection_density_floor)
                    beta_face=this%time_step/rho_face
                    radius_face=mesh%mesh(1,loop(1,1),j,k)-0.5_dp*dx(1)+real(ii-1,dp)*dx(1)
                    metric_face=1.0_dp
                    if(geometry_power>0)metric_face=max(radius_face,0.0_dp)**geometry_power
                    conductance_x(ii,jj,kk)=beta_face*base_volume*metric_face/dx(1)**2
                    il=i-1
                    jl=j
                    kl=k
                    ir=i
                    jr=j
                    kr=k
                    bound_l=bc%bc_markers(il,jl,kl)
                    bound_r=bc%bc_markers(ir,jr,kr)
                    if(bound_l/=0.or.bound_r/=0)then
                        boundary_marker=merge(bound_l,bound_r,bound_l/=0)
                        boundary_name=bc%boundary_types(boundary_marker)%get_type_name()
                        elliptic_boundary%type_x(ii,jj,kk)=merge(elliptic_bc_dirichlet,elliptic_bc_neumann, &
                            boundary_name=='outlet')
                    end if
                end do
            end do
        end do

        if(dimensions>=2)then
            do kk=1,size(solution,3)
                do jj=1,size(solution,2)+1
                    do ii=1,size(solution,1)
                        i=loop(1,1)+ii-1
                        j=loop(2,1)+jj-1
                        k=loop(3,1)+kk-1
                        rho_face=max(density_face_old%cells(2,i,j,k),cabaret_projection_density_floor)
                        beta_face=this%time_step/rho_face
                        radius=mesh%mesh(1,i,loop(2,1),k)
                        metric_center=1.0_dp
                        if(geometry_power>0)metric_center=max(radius,0.0_dp)**geometry_power
                        conductance_y(ii,jj,kk)=beta_face*base_volume*metric_center/dx(2)**2
                        il=i
                        jl=j-1
                        kl=k
                        ir=i
                        jr=j
                        kr=k
                        bound_l=bc%bc_markers(il,jl,kl)
                        bound_r=bc%bc_markers(ir,jr,kr)
                        if(bound_l/=0.or.bound_r/=0)then
                            boundary_marker=merge(bound_l,bound_r,bound_l/=0)
                            boundary_name=bc%boundary_types(boundary_marker)%get_type_name()
                            elliptic_boundary%type_y(ii,jj,kk)=merge(elliptic_bc_dirichlet,elliptic_bc_neumann, &
                                boundary_name=='outlet')
                        end if
                    end do
                end do
            end do
        end if
        if(dimensions>=3)then
            do kk=1,size(solution,3)+1
                do jj=1,size(solution,2)
                    do ii=1,size(solution,1)
                        i=loop(1,1)+ii-1
                        j=loop(2,1)+jj-1
                        k=loop(3,1)+kk-1
                        rho_face=max(density_face_old%cells(3,i,j,k),cabaret_projection_density_floor)
                        beta_face=this%time_step/rho_face
                        radius=mesh%mesh(1,i,j,loop(3,1))
                        metric_center=1.0_dp
                        if(geometry_power>0)metric_center=max(radius,0.0_dp)**geometry_power
                        conductance_z(ii,jj,kk)=beta_face*base_volume*metric_center/dx(3)**2
                        il=i
                        jl=j
                        kl=k-1
                        ir=i
                        jr=j
                        kr=k
                        bound_l=bc%bc_markers(il,jl,kl)
                        bound_r=bc%bc_markers(ir,jr,kr)
                        if(bound_l/=0.or.bound_r/=0)then
                            boundary_marker=merge(bound_l,bound_r,bound_l/=0)
                            boundary_name=bc%boundary_types(boundary_marker)%get_type_name()
                            elliptic_boundary%type_z(ii,jj,kk)=merge(elliptic_bc_dirichlet,elliptic_bc_neumann, &
                                boundary_name=='outlet')
                        end if
                    end do
                end do
            end do
        end if
    end associate

    this%pressure_solver%relative_tolerance = cabaret_projection_relative_tolerance
    this%pressure_solver%absolute_tolerance=cabaret_projection_absolute_tolerance
    call this%pressure_solver%solve(solution,rhs,conductance_x,conductance_y,conductance_z,cell_volume,active, &
        elliptic_boundary,dimensions,.false.,statistics)
    if(.not.statistics%converged)error stop 'CABARET multigrid projection did not converge'

    associate(delta => this%pressure_correction%s_ptr, bc => this%boundary%bc_ptr)
        delta%cells=0.0_dp
        do k=loop(3,1),loop(3,2)
            do j=loop(2,1),loop(2,2)
                do i=loop(1,1),loop(1,2)
                    if(bc%bc_markers(i,j,k)/=0)cycle
                    ii=i-loop(1,1)+1
                    jj=j-loop(2,1)+1
                    kk=k-loop(3,1)+1
                    delta%cells(i,j,k)=solution(ii,jj,kk)
                    if(accumulate_pressure)then
                        this%p_dyn%s_ptr%cells(i,j,k)=this%p_dyn%s_ptr%cells(i,j,k)+solution(ii,jj,kk)
                    else
                        this%p_dyn%s_ptr%cells(i,j,k)=solution(ii,jj,kk)
                    end if
                end do
            end do
        end do
        call this%mpi_support%exchange_conservative_scalar_field(delta)
        call this%mpi_support%exchange_conservative_scalar_field(this%p_dyn%s_ptr)

        do dim=1,dimensions
            face_loop=flow_loop
            do dim1=1,dimensions
                face_loop(dim1,2)=flow_loop(dim1,2)-(1-I_m(dim1,dim))
            end do
            do k=face_loop(3,1),face_loop(3,2)
                do j=face_loop(2,1),face_loop(2,2)
                    do i=face_loop(1,1),face_loop(1,2)
                        il=i-I_m(dim,1)
                        jl=j-I_m(dim,2)
                        kl=k-I_m(dim,3)
                        ir=i
                        jr=j
                        kr=k
                        bound_l=bc%bc_markers(il,jl,kl)
                        bound_r=bc%bc_markers(ir,jr,kr)
                        if(bound_l==0.and.bound_r==0)then
                            grad_phi=(delta%cells(ir,jr,kr)-delta%cells(il,jl,kl))/dx(dim)
                        else if(bound_l/=0.and.bound_r==0)then
                            type_l=bc%boundary_types(bound_l)%get_type_name()
                            if(type_l/='outlet')cycle
                            grad_phi=2.0_dp*delta%cells(ir,jr,kr)/dx(dim)
                        else if(bound_l==0.and.bound_r/=0)then
                            type_r=bc%boundary_types(bound_r)%get_type_name()
                            if(type_r/='outlet')cycle
                            grad_phi=-2.0_dp*delta%cells(il,jl,kl)/dx(dim)
                        else
                            cycle
                        end if
                        rho_face=max(this%rho_f_new%s_ptr%cells(dim,i,j,k),cabaret_projection_density_floor)
                        this%v_f_new%v_ptr%pr(dim)%cells(dim,i,j,k)=this%v_f_new%v_ptr%pr(dim)%cells(dim,i,j,k)- &
                            this%time_step*grad_phi/rho_face
                    end do
                end do
            end do
        end do
    end associate
    call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
    call this%update_projection_pressure_faces()
end subroutine apply_low_mach_pressure_projection

subroutine update_projection_pressure_faces(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: dimensions,dim,dim1,i,j,k,il,jl,kl,ir,jr,kr,bound_l,bound_r
    integer, dimension(3,2) :: flow_loop,face_loop
    real(dp) :: phi_face
    character(len=20) :: type_l,type_r
    dimensions=this%domain%get_domain_dimensions()
    flow_loop=this%domain%get_local_inner_faces_bounds()
    associate(phi=>this%p_dyn%s_ptr,bc=>this%boundary%bc_ptr,pface=>this%p_f_new%s_ptr)
        do dim=1,dimensions
            face_loop=flow_loop
            do dim1=1,dimensions
                face_loop(dim1,2)=flow_loop(dim1,2)-(1-I_m(dim1,dim))
            end do
            do k=face_loop(3,1),face_loop(3,2)
                do j=face_loop(2,1),face_loop(2,2)
                    do i=face_loop(1,1),face_loop(1,2)
                        il=i-I_m(dim,1)
                        jl=j-I_m(dim,2)
                        kl=k-I_m(dim,3)
                        ir=i
                        jr=j
                        kr=k
                        bound_l=bc%bc_markers(il,jl,kl)
                        bound_r=bc%bc_markers(ir,jr,kr)
                        if (bound_l==0 .and. bound_r==0) then
                            phi_face=0.5_dp*(phi%cells(il,jl,kl)+phi%cells(ir,jr,kr))
                        else if (bound_l/=0 .and. bound_r==0) then
                            type_l=bc%boundary_types(bound_l)%get_type_name()
                            phi_face=merge(0.0_dp,phi%cells(ir,jr,kr),type_l=='outlet')
                        else if (bound_l==0 .and. bound_r/=0) then
                            type_r=bc%boundary_types(bound_r)%get_type_name()
                            phi_face=merge(0.0_dp,phi%cells(il,jl,kl),type_r=='outlet')
                        else
                            phi_face=0.0_dp
                        end if
                        this%dynamic_pressure_face(dim,i,j,k)=phi_face
                        pface%cells(dim,i,j,k)=this%thermodynamic_pressure+phi_face
                    end do
                end do
            end do
        end do
    end associate
    call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
end subroutine update_projection_pressure_faces

subroutine correct_projection_cell_velocity(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: dimensions,i,j,k,dim
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: grad_phi,rho_cell
    dimensions=this%domain%get_domain_dimensions()
    loop=this%domain%get_local_inner_cells_bounds()
    dx=this%mesh%mesh_ptr%get_cell_edges_length()
    associate(v => this%v%v_ptr, rho => this%rho%s_ptr, bc => this%boundary%bc_ptr)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    rho_cell = max(rho%cells(i,j,k), cabaret_projection_density_floor)
                    do dim = 1, dimensions
                        grad_phi = (this%dynamic_pressure_face(dim, i+I_m(dim,1), j+I_m(dim,2), k+I_m(dim,3)) - &
                            this%dynamic_pressure_face(dim,i,j,k))/dx(dim)
                        v%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k) - &
                            this%time_step*grad_phi/rho_cell
                    end do
                end do
            end do
        end do
    end associate
    call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)
end subroutine correct_projection_cell_velocity

! Evaluate continuum physics once at the old state and assemble rates for
! momentum, sensible enthalpy, species, and the divergence constraint.
subroutine evaluate_continuum_source_rates(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, species_number
    integer :: i, j, k, dim, spec, phase
    real(dp) :: particle_mass_rate, particle_energy_rate
    real(dp) :: particle_kinetic_energy, particle_momentum_rate
    integer, dimension(3,2) :: loop

    dimensions = this%domain%get_domain_dimensions()
    species_number = this%chem%chem_ptr%species_number
    loop = this%domain%get_local_inner_cells_bounds()

    this%momentum_source = 0.0_dp
    this%sensible_enthalpy_source = 0.0_dp
    this%species_mass_source = 0.0_dp

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

    call this%advance_particle_phases()

    associate(bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(guided) private(i,j,k,dim,spec,phase,particle_mass_rate,particle_energy_rate, &
        !$omp& particle_kinetic_energy,particle_momentum_rate)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle

                    if (this%reactive_flag) then
                        this%sensible_enthalpy_source(i,j,k) = this%sensible_enthalpy_source(i,j,k)+ &
                            this%E_f_prod_chem%s_ptr%cells(i,j,k)
                        do spec = 1, species_number
                            this%species_mass_source(spec,i,j,k) = this%species_mass_source(spec,i,j,k)+ &
                                this%Y_prod_chem%v_ptr%pr(spec)%cells(i,j,k)
                        end do
                    end if

                    if (this%heat_trans_flag) then
                        this%sensible_enthalpy_source(i,j,k) = this%sensible_enthalpy_source(i,j,k)+ &
                            this%E_f_prod_heat%s_ptr%cells(i,j,k)
                    end if

                    if (this%radiation_flag) then
                        this%sensible_enthalpy_source(i,j,k) = this%sensible_enthalpy_source(i,j,k)+ &
                            this%E_f_prod_rad%s_ptr%cells(i,j,k)
                    end if

                    if (this%diffusion_flag) then
                        this%sensible_enthalpy_source(i,j,k) = this%sensible_enthalpy_source(i,j,k)+ &
                            this%E_f_prod_diff%s_ptr%cells(i,j,k)
                        do spec = 1, species_number
                            this%species_mass_source(spec,i,j,k) = this%species_mass_source(spec,i,j,k)+ &
                                this%Y_prod_diff%v_ptr%pr(spec)%cells(i,j,k)
                        end do
                    end if

                    if (this%viscosity_flag) then
                        do dim = 1, dimensions
                            this%momentum_source(dim,i,j,k) = this%momentum_source(dim,i,j,k)+ &
                                this%v_prod_visc%v_ptr%pr(dim)%cells(i,j,k)
                        end do
                    end if

                    ! Reduced-gravity buoyancy in the FDS velocity convention.  FDS stores
                    ! F_a=(rho-buoyancy_reference_density)g/rho and updates u <- u-dt*F_a, so the actual
                    ! acceleration is (buoyancy_reference_density-rho)g/rho.
                    do dim = 1, dimensions
                        this%momentum_source(dim,i,j,k) = this%momentum_source(dim,i,j,k)+ &
                            (this%buoyancy_reference_density-this%rho_old(i,j,k))*this%gravity(dim)
                    end do

                    ! Convert the common dispersed-phase source contract to the
                    ! low-Mach conservative momentum and sensible-enthalpy rates.
                    do phase = 1, this%additional_particles_phases_number
                        particle_mass_rate = &
                            this%rho_prod_particles(phase)%s_ptr%cells(i,j,k)
                        particle_energy_rate = &
                            this%E_f_prod_particles(phase)%s_ptr%cells(i,j,k)
                        particle_kinetic_energy = 0.0_dp
                        do dim = 1, dimensions
                            particle_momentum_rate = this%rho_old(i,j,k) * &
                                this%v_prod_particles(phase)%v_ptr%pr(dim)%cells(i,j,k)
                            this%momentum_source(dim,i,j,k) = &
                                this%momentum_source(dim,i,j,k) + particle_momentum_rate
                            particle_energy_rate = particle_energy_rate - &
                                this%v_old(dim,i,j,k)*particle_momentum_rate
                            particle_kinetic_energy = particle_kinetic_energy + &
                                0.5_dp*this%v_old(dim,i,j,k)**2
                        end do
                        this%sensible_enthalpy_source(i,j,k) = &
                            this%sensible_enthalpy_source(i,j,k) + &
                            particle_energy_rate + &
                            particle_kinetic_energy*particle_mass_rate
                        do spec = 1, species_number
                            this%species_mass_source(spec,i,j,k) = &
                                this%species_mass_source(spec,i,j,k) + &
                                this%Y_prod_particles(phase)%v_ptr%pr(spec)%cells(i,j,k)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end associate
end subroutine evaluate_continuum_source_rates

	subroutine material_characteristic_source_rates(this,rho_state,T_state,Y_state,Q_h,Q_Y,g_Y,g_T)
		class(cabaret_low_mach_solver),intent(in)::this
		real(dp),intent(in)::rho_state,T_state,Q_h
		real(dp),dimension(:),intent(in)::Y_state,Q_Y
		real(dp),dimension(:),intent(out)::g_Y
		real(dp),intent(out)::g_T
		integer::spec
		real(dp)::M,cp_molar,cp_mass,sum_hq,total_mass_source
		real(dp),dimension(size(Y_state))::X,hsk

		if (rho_state<=tiny(1.0_dp)) &
			error stop 'CABARET material source shift: non-positive density'

		M=this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_state)
		total_mass_source=sum(Q_Y)
		do spec=1,size(Y_state)
			X(spec)=Y_state(spec)*M/this%thermo%thermo_ptr%molar_masses(spec)
			hsk(spec)=(this%thermo%thermo_ptr%specie_enthalpy_molar(T_state,spec)- &
				this%thermo%thermo_ptr%specie_enthalpy_molar(T_ref,spec))/ &
				this%thermo%thermo_ptr%molar_masses(spec)
			! Exact mass-fraction source for non-zero dispersed-phase mass transfer.
			g_Y(spec)=(Q_Y(spec)-Y_state(spec)*total_mass_source)/rho_state
		end do

		cp_molar=this%thermo%thermo_ptr%mixture_cp_molar(T_state,X)
		cp_mass=cp_molar/M
		if (cp_mass<=tiny(1.0_dp)) &
			error stop 'CABARET material source shift: non-positive heat capacity'

		sum_hq=0.0_dp
		do spec=1,size(Y_state)
			sum_hq=sum_hq+hsk(spec)*Q_Y(spec)
		end do
		g_T=(Q_h-sum_hq)/(rho_state*cp_mass)
	end subroutine material_characteristic_source_rates

subroutine update_cell_thermodynamics(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    call cabaret_eos_timer%tic()
    call this%state_eq%apply_state_equation_low_mach_rsst_enthalpy(this%thermodynamic_pressure)
    call cabaret_eos_timer%toc(new_iter=.true.)
end subroutine update_cell_thermodynamics

subroutine set_face_thermodynamic_pressure(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, dim, dim1, i, j, k, il, jl, kl
    integer, dimension(3,2) :: flow_inner_loop, loop

    dimensions = this%domain%get_domain_dimensions()
    flow_inner_loop = this%domain%get_local_inner_faces_bounds()

    associate(density_face_old => this%rho_f_new%s_ptr, p_f => this%p_f_new%s_ptr, &
              bc => this%boundary%bc_ptr)
    do dim = 1, dimensions
        loop = flow_inner_loop
        do dim1 = 1, dimensions
            loop(dim1,2) = flow_inner_loop(dim1,2) - (1-I_m(dim1,dim))
        end do
        !$omp parallel do collapse(3) schedule(static) private(i,j,k,il,jl,kl)
        do k = loop(3,1), loop(3,2)
        do j = loop(2,1), loop(2,2)
        do i = loop(1,1), loop(1,2)
            il = i-I_m(dim,1)
            jl = j-I_m(dim,2)
            kl = k-I_m(dim,3)
            if ((bc%bc_markers(i,j,k) == 0) .or. (bc%bc_markers(il,jl,kl) == 0)) then
                if (density_face_old%cells(dim,i,j,k) <= 0.0_dp) &
                    error stop 'CABARET projection: non-positive face density'
                p_f%cells(dim,i,j,k) = this%thermodynamic_pressure
            end if
        end do
        end do
        end do
        !$omp end parallel do
    end do
    end associate

    call this%exchange_face_pressure_density()
end subroutine set_face_thermodynamic_pressure

subroutine update_flow_thermodynamics(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    call cabaret_eos_timer%tic()
    call this%set_face_thermodynamic_pressure()
    call this%state_eq%apply_state_equation_flow_variables()
    this%p_f_new%s_ptr%cells = this%thermodynamic_pressure + this%dynamic_pressure_face
    call this%exchange_face_pressure_density()
    call cabaret_eos_timer%toc(new_iter=.true.)
end subroutine update_flow_thermodynamics

subroutine finalize_material_face_reconstruction(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    call this%normalize_face_mass_fractions()
    call this%exchange_face_primitive_state()
end subroutine finalize_material_face_reconstruction

	real(dp) function temperature_face_value(this,rho_face,Y_face) result(T_face)
		class(cabaret_low_mach_solver),intent(in)::this
		real(dp),intent(in)::rho_face
		real(dp),dimension(:),intent(in)::Y_face
		real(dp)::M_face

		if (rho_face<=0.0_dp) &
			error stop 'CABARET temperature reconstruction: non-positive face density'
		M_face=this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_face)
		T_face=this%thermodynamic_pressure*M_face/(rho_face*r_gase_J)
		if (T_face<=1.0_dp .or. T_face>=tables_temperature_ceiling) &
			error stop 'CABARET temperature reconstruction: old face temperature outside table'
	end function temperature_face_value

	real(dp) function sensible_enthalpy_face_value(this,rho_face,Y_face) result(h_face)
		class(cabaret_low_mach_solver),intent(in)::this
		real(dp),intent(in)::rho_face
		real(dp),dimension(:),intent(in)::Y_face
		integer::spec
		real(dp)::T_face,M_face,h_molar,h_ref_molar
		real(dp),dimension(size(Y_face))::X_face

		T_face=this%temperature_face_value(rho_face,Y_face)
		M_face=this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Y_face)
		do spec=1,size(Y_face)
			X_face(spec)=Y_face(spec)*M_face/this%thermo%thermo_ptr%molar_masses(spec)
		end do
		h_molar=this%thermo%thermo_ptr%mixture_enthalpy_molar(T_face,X_face)
		h_ref_molar=this%thermo%thermo_ptr%mixture_enthalpy_molar(T_ref,X_face)
		h_face=(h_molar-h_ref_molar)/M_face
	end function sensible_enthalpy_face_value

	! Close the reconstructed contact state with the low-Mach thermodynamic
	! pressure.  T_f and Y_f are independent CABARET material variables;
	! density_face_old is the dependent EOS quantity required by the mass flux.
	subroutine close_material_face_density(this)
		class(cabaret_low_mach_solver),intent(inout)::this
		integer::dimensions,species_number,dim,dim1,spec,i,j,k,il,jl,kl
		integer,dimension(3,2)::flow_inner_loop,loop
		real(dp) :: M_face, T_face
		real(dp), dimension(:), allocatable :: Y_face

		dimensions=this%domain%get_domain_dimensions()
		species_number=this%chem%chem_ptr%species_number
		flow_inner_loop=this%domain%get_local_inner_faces_bounds()

		associate(rho_f_new => this%rho_f_new%s_ptr, &
			T_f_new => this%T_f_new%s_ptr, Y_f_new => this%Y_f_new%v_ptr, &
			bc => this%boundary%bc_ptr)
		!$omp parallel default(shared) private(dim,dim1,spec,i,j,k,il,jl,kl,loop, &
		!$omp& M_face,T_face,Y_face)
		allocate(Y_face(species_number))
		do dim=1,dimensions
			loop=flow_inner_loop
			do dim1=1,dimensions
				loop(dim1,2)=flow_inner_loop(dim1,2)-(1-I_m(dim1,dim))
			end do
			!$omp do collapse(3) schedule(guided)
			do k=loop(3,1),loop(3,2)
			do j=loop(2,1),loop(2,2)
			do i=loop(1,1),loop(1,2)
				il=i-I_m(dim,1);jl=j-I_m(dim,2);kl=k-I_m(dim,3)
				if ((bc%bc_markers(i,j,k)/=0) .and. &
					(bc%bc_markers(il,jl,kl)/=0)) cycle

				T_face=T_f_new%cells(dim,i,j,k)
				if (T_face<=1.0_dp .or. T_face>=tables_temperature_ceiling) then
					write(*,*) 'CABARET temperature face closure failure ', &
						dim,i,j,k,T_face
					error stop 'CABARET temperature face closure: T outside table'
				end if

				do spec=1,species_number
					Y_face(spec)=Y_f_new%pr(spec)%cells(dim,i,j,k)
				end do
				M_face=this%thermo%thermo_ptr% &
					mixture_molar_mass_from_mass_fractions(Y_face)
				rho_f_new%cells(dim,i,j,k) = this%thermodynamic_pressure*M_face/(r_gase_J*T_face)
			end do
			end do
			end do
			!$omp end do nowait
		end do
		deallocate(Y_face)
		!$omp end parallel
		end associate

		call this%exchange_face_pressure_density()
	end subroutine close_material_face_density

subroutine predict_material_half_step(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: dimensions, species_number, nu
    integer :: i,j,k,dim,spec
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: r, w_l, w_c, w_r, u_l, u_r
    real(dp) :: rho_l, rho_r, rho_c, flux_l, flux_r
    real(dp) :: hs_l, hs_r, rhs_h, sum_y, rhs_y, pressure_work_rate, total_mass_source
    real(dp), dimension(:), allocatable :: Y_l, Y_r, Y_new

    dimensions=this%domain%get_domain_dimensions()
    species_number=this%chem%chem_ptr%species_number
    loop=this%domain%get_local_inner_cells_bounds()
    dx=this%mesh%mesh_ptr%get_cell_edges_length()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        nu=1
    case('cylindrical')
        nu=2
    case('spherical')
        nu=3
    case default
        nu=1
    end select

    associate(hs => this%h_s%s_ptr, Y => this%Y%v_ptr, &
        velocity_face_old => this%velocity_face_old, &
        species_face_old => this%species_face_old, &
        density_face_old => this%density_face_old, &
        bc => this%boundary%bc_ptr)
        !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,spec,r,w_l,w_c,w_r,u_l,u_r,rho_l,rho_r,rho_c,flux_l,flux_r, &
            !$omp& hs_l,hs_r,rhs_h,sum_y,rhs_y,pressure_work_rate,total_mass_source,Y_l,Y_r,Y_new)
        allocate(Y_l(species_number),Y_r(species_number),Y_new(species_number))
        !$omp do collapse(3) schedule(guided)
        do k=loop(3,1),loop(3,2)
            do j=loop(2,1),loop(2,2)
                do i=loop(1,1),loop(1,2)
                    if (bc%bc_markers(i,j,k)/=0) cycle
                    rhs_h=0.0_dp
                    Y_new=0.0_dp
                    r= this%mesh%mesh_ptr%mesh(1,i,j,k)
                    rho_c=this%rho_old(i,j,k)
                    if (rho_c<=tiny(1.0_dp)) error stop 'CABARET predictor: non-positive cell density'
                    do dim=1,dimensions
                        w_l=1.0_dp
                        w_c=1.0_dp
                        w_r=1.0_dp
                        if (dim==1 .and. nu>1) then
                            w_l=(r-0.5_dp*dx(1))**(nu-1)
                            w_c=r**(nu-1)
                            w_r=(r+0.5_dp*dx(1))**(nu-1)
                        end if
                        u_l=velocity_face_old(dim,dim,i,j,k)
                        u_r=velocity_face_old(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        rho_l=density_face_old(dim,i,j,k)
                        rho_r=density_face_old(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        if ((rho_l<=tiny(1.0_dp)) .or. (rho_r<=tiny(1.0_dp))) &
                            error stop 'CABARET predictor: non-positive face density'
                        flux_l=rho_l*u_l
                        flux_r=rho_r*u_r
                        Y_l(:)=species_face_old(:,dim,i,j,k)
                        Y_r(:)=species_face_old(:,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        ! The old face state has no separate persistent h_s array.  Recovering
                        ! h_s from its p0-consistent rho,Y state is algebraically identical and
                        ! also initializes the very first physical step safely.
                        hs_l=this%sensible_enthalpy_face_value(rho_l,Y_l)
                        hs_r=this%sensible_enthalpy_face_value(rho_r,Y_r)
                        rhs_h=rhs_h-(flux_r*w_r*(hs_r-this%h_s_old(i,j,k))- &
                            flux_l*w_l*(hs_l-this%h_s_old(i,j,k)))/(rho_c*dx(dim)*w_c)
                        do spec=1,species_number
                            Y_new(spec)=Y_new(spec)-(flux_r*w_r*(Y_r(spec)-this%Y_old(spec,i,j,k))- &
                                flux_l*w_l*(Y_l(spec)-this%Y_old(spec,i,j,k)))/(rho_c*dx(dim)*w_c)
                        end do
                    end do
                    this%sensible_enthalpy_advection_old(i,j,k)=rhs_h
                    pressure_work_rate = 0.0_dp
                    if (this%closed_domain) &
                        pressure_work_rate = this%thermodynamic_pressure_rate_old
                    hs%cells(i,j,k)=this%h_s_old(i,j,k)+0.5_dp*this%time_step* &
                        (rhs_h+(this%sensible_enthalpy_source(i,j,k)+pressure_work_rate)/this%rho_old(i,j,k))
                    total_mass_source=sum(this%species_mass_source(:,i,j,k))
                    sum_y=0.0_dp
                    do spec=1,species_number
                        this%species_advection_old(spec,i,j,k)=Y_new(spec)
                        rhs_y=Y_new(spec)+(this%species_mass_source(spec,i,j,k)- &
                            this%Y_old(spec,i,j,k)*total_mass_source)/this%rho_old(i,j,k)
                        Y_new(spec)=max(this%Y_old(spec,i,j,k)+0.5_dp*this%time_step*rhs_y,0.0_dp)
                        sum_y=sum_y+Y_new(spec)
                    end do
                    if (sum_y<=tiny(1.0_dp)) error stop 'CABARET predictor: zero species sum'
                    do spec=1,species_number
                        Y%pr(spec)%cells(i,j,k)=Y_new(spec)/sum_y
                    end do
                end do
            end do
        end do
        !$omp end do nowait
        deallocate(Y_l,Y_r,Y_new)
        !$omp end parallel
    end associate
end subroutine predict_material_half_step

subroutine correct_material_full_step(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: dimensions, species_number, nu
    integer :: i,j,k,dim,spec
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: r,w_l,w_c,w_r,u_l,u_r,rho_l,rho_r,rho_c,flux_l,flux_r
    real(dp) :: hs_l,hs_r,rhs_h,sum_y,rhs_y,Y_raw,pressure_work_rate,total_mass_source
    real(dp), dimension(:), allocatable :: Y_l,Y_r,Y_new

    dimensions=this%domain%get_domain_dimensions()
    species_number=this%chem%chem_ptr%species_number
    loop=this%domain%get_local_inner_cells_bounds()
    dx=this%mesh%mesh_ptr%get_cell_edges_length()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        nu=1
    case('cylindrical')
        nu=2
    case('spherical')
        nu=3
    case default
        nu=1
    end select
    call this%exchange_face_thermodynamic_state()
    associate(hs=>this%h_s%s_ptr,Y=>this%Y%v_ptr,rho=>this%rho%s_ptr,density_face_old=>this%rho_f_new%s_ptr, &
            velocity_face_old=>this%v_f_new%v_ptr,species_face_old=>this%Y_f_new%v_ptr,bc=>this%boundary%bc_ptr)
        !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,spec,r,w_l,w_c,w_r,u_l,u_r,rho_l,rho_r,rho_c,flux_l,flux_r, &
            !$omp& hs_l,hs_r,rhs_h,sum_y,rhs_y,Y_raw,pressure_work_rate,total_mass_source,Y_l,Y_r,Y_new)
        allocate(Y_l(species_number),Y_r(species_number),Y_new(species_number))
        !$omp do collapse(3) schedule(guided)
        do k=loop(3,1),loop(3,2)
            do j=loop(2,1),loop(2,2)
                do i=loop(1,1),loop(1,2)
                    if (bc%bc_markers(i,j,k)/=0) cycle
                    rhs_h=0.0_dp
                    Y_new=0.0_dp
                    r=this%mesh%mesh_ptr%mesh(1,i,j,k)
                    rho_c=rho%cells(i,j,k)
                    if (rho_c<=tiny(1.0_dp)) error stop 'CABARET corrector: non-positive cell density'
                    do dim=1,dimensions
                        w_l=1.0_dp
                        w_c=1.0_dp
                        w_r=1.0_dp
                        if (dim==1 .and. nu>1) then
                            w_l=(r-0.5_dp*dx(1))**(nu-1)
                            w_c=r**(nu-1)
                            w_r=(r+0.5_dp*dx(1))**(nu-1)
                        end if
                        u_l=velocity_face_old%pr(dim)%cells(dim,i,j,k)
                        u_r=velocity_face_old%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        rho_l=density_face_old%cells(dim,i,j,k)
                        rho_r=density_face_old%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        if ((rho_l<=tiny(1.0_dp)) .or. (rho_r<=tiny(1.0_dp))) &
                            error stop 'CABARET corrector: non-positive face density'
                        flux_l=rho_l*u_l
                        flux_r=rho_r*u_r
                        do spec=1,species_number
                            Y_l(spec)=species_face_old%pr(spec)%cells(dim,i,j,k)
                            Y_r(spec)=species_face_old%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
                        end do
                        ! density_face_old,species_face_old have already been closed from the limited h_s invariant
                        ! and exchanged across MPI boundaries.  Recovering h_s here is therefore
                        ! equivalent to the stored invariant and avoids a separate face-h_s halo.
                        hs_l=this%sensible_enthalpy_face_value(rho_l,Y_l)
                        hs_r=this%sensible_enthalpy_face_value(rho_r,Y_r)
                        rhs_h=rhs_h-(flux_r*w_r*(hs_r-hs%cells(i,j,k))- &
                            flux_l*w_l*(hs_l-hs%cells(i,j,k)))/(rho_c*dx(dim)*w_c)
                        do spec=1,species_number
                            Y_new(spec)=Y_new(spec)-(flux_r*w_r*(Y_r(spec)-Y%pr(spec)%cells(i,j,k))- &
                                flux_l*w_l*(Y_l(spec)-Y%pr(spec)%cells(i,j,k)))/(rho_c*dx(dim)*w_c)
                        end do
                    end do
                    pressure_work_rate = 0.0_dp
                    if (this%closed_domain) &
                        pressure_work_rate = this%thermodynamic_pressure_rate_mid
                    hs%cells(i,j,k)=this%h_s_old(i,j,k)+0.5_dp*this%time_step* &
                        (this%sensible_enthalpy_advection_old(i,j,k)+rhs_h)+this%time_step* &
                        (this%sensible_enthalpy_source(i,j,k)+pressure_work_rate)/this%rho_old(i,j,k)
                    total_mass_source=sum(this%species_mass_source(:,i,j,k))
                    sum_y=0.0_dp
                    do spec=1,species_number
                        rhs_y=(this%species_mass_source(spec,i,j,k)- &
                            this%Y_old(spec,i,j,k)*total_mass_source)/this%rho_old(i,j,k)
                        Y_raw=this%Y_old(spec,i,j,k)+0.5_dp*this%time_step* &
                            (this%species_advection_old(spec,i,j,k)+Y_new(spec))+this%time_step*rhs_y
                        Y_new(spec)=max(Y_raw,0.0_dp)
                        sum_y=sum_y+Y_new(spec)
                    end do
                    if (sum_y<=tiny(1.0_dp)) error stop 'CABARET corrector: zero species sum'
                    do spec=1,species_number
                        Y_raw=Y_new(spec)/sum_y
                        Y%pr(spec)%cells(i,j,k)=Y_raw
                    end do
                end do
            end do
        end do
        !$omp end do nowait
        deallocate(Y_l,Y_r,Y_new)
        !$omp end parallel
    end associate
end subroutine correct_material_full_step

subroutine predict_momentum_half_step(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, geometry_power
    integer :: i, j, k, dim, direction, ih, jh, kh
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: radius, flux_high, flux_low

    dimensions = this%domain%get_domain_dimensions()
    loop = this%domain%get_local_inner_cells_bounds()
    dx = this%mesh%mesh_ptr%get_cell_edges_length()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        geometry_power = 1
    case('cylindrical')
        geometry_power = 2
    case('spherical')
        geometry_power = 3
    case default
        error stop 'CABARET projection: unsupported coordinate system'
    end select

    associate(rho => this%rho%s_ptr,v => this%v%v_ptr,bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(guided) &
            !$omp& private(i,j,k,dim,direction,ih,jh,kh,radius,flux_high,flux_low)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
                    do dim = 1, dimensions
                        v%pr(dim)%cells(i,j,k) = 0.0_dp
                        do direction = 1, dimensions
                            ih = i+I_m(direction,1)
                            jh = j+I_m(direction,2)
                            kh = k+I_m(direction,3)
                            flux_high = this%density_face_old(direction,ih,jh,kh)* &
                                this%velocity_face_old(dim,direction,ih,jh,kh)*this%velocity_face_old(direction,direction,ih,jh,kh)
                            flux_low = this%density_face_old(direction,i,j,k)* &
                                this%velocity_face_old(dim,direction,i,j,k)*this%velocity_face_old(direction,direction,i,j,k)
                            v%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k)- &
                                (flux_high-flux_low)/dx(direction)
                            if (direction == 1 .and. geometry_power > 1) &
                                v%pr(dim)%cells(i,j,k) = v%pr(dim)%cells(i,j,k)- &
                                finite_volume_geometry_coefficient(geometry_power,radius,dx(1))* &
                                0.5_dp*(flux_high+flux_low)
                        end do
                        this%momentum_advection_old(dim,i,j,k) = v%pr(dim)%cells(i,j,k)
                        v%pr(dim)%cells(i,j,k) = (this%rho_old(i,j,k)*this%v_old(dim,i,j,k)+ &
                            0.5_dp*this%time_step*(this%momentum_advection_old(dim,i,j,k)+ &
                            this%momentum_source(dim,i,j,k)))/rho%cells(i,j,k)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end associate
end subroutine predict_momentum_half_step

subroutine correct_momentum_full_step(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, geometry_power
    integer :: i, j, k, dim, direction, ih, jh, kh
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: radius, flux_high, flux_low, rhs_new

    dimensions = this%domain%get_domain_dimensions()
    loop = this%domain%get_local_inner_cells_bounds()
    dx = this%mesh%mesh_ptr%get_cell_edges_length()
    select case(this%domain%get_coordinate_system_name())
    case('cartesian')
        geometry_power = 1
    case('cylindrical')
        geometry_power = 2
    case('spherical')
        geometry_power = 3
    case default
        error stop 'CABARET projection: unsupported coordinate system'
    end select

    call this%exchange_face_thermodynamic_state()
    associate(rho => this%rho%s_ptr,v => this%v%v_ptr,bc => this%boundary%bc_ptr)
        !$omp parallel do collapse(3) schedule(guided) &
            !$omp& private(i,j,k,dim,direction,ih,jh,kh,radius,flux_high,flux_low,rhs_new)
        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    radius = this%mesh%mesh_ptr%mesh(1,i,j,k)
                    do dim = 1, dimensions
                        rhs_new = 0.0_dp
                        do direction = 1, dimensions
                            ih = i+I_m(direction,1)
                            jh = j+I_m(direction,2)
                            kh = k+I_m(direction,3)
                            flux_high = this%rho_f_new%s_ptr%cells(direction,ih,jh,kh)* &
                                this%v_f_new%v_ptr%pr(dim)%cells(direction,ih,jh,kh)* &
                                this%v_f_new%v_ptr%pr(direction)%cells(direction,ih,jh,kh)
                            flux_low = this%rho_f_new%s_ptr%cells(direction,i,j,k)* &
                                this%v_f_new%v_ptr%pr(dim)%cells(direction,i,j,k)* &
                                this%v_f_new%v_ptr%pr(direction)%cells(direction,i,j,k)
                            rhs_new = rhs_new-(flux_high-flux_low)/dx(direction)
                            if (direction == 1 .and. geometry_power > 1) &
                                rhs_new = rhs_new- &
                                finite_volume_geometry_coefficient(geometry_power,radius,dx(1))* &
                                0.5_dp*(flux_high+flux_low)
                        end do
                        v%pr(dim)%cells(i,j,k) = (this%rho_old(i,j,k)*this%v_old(dim,i,j,k)+ &
                            0.5_dp*this%time_step*(this%momentum_advection_old(dim,i,j,k)+rhs_new)+ &
                            this%time_step*this%momentum_source(dim,i,j,k))/rho%cells(i,j,k)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end associate
end subroutine correct_momentum_full_step

subroutine refresh_thermodynamic_ghosts(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: i,j,k,spec,species_number
    integer,dimension(3,2)::loop
    real(dp)::M,h,h0,sumy
    real(dp),dimension(:),allocatable::X,Yloc
    species_number=this%chem%chem_ptr%species_number
    loop=this%domain%get_local_utter_cells_bounds()
    allocate(X(species_number),Yloc(species_number))
    associate(T=>this%T%s_ptr,rho=>this%rho%s_ptr,Y=>this%Y%v_ptr,hs=>this%h_s%s_ptr,bc=>this%boundary%bc_ptr)
        do k=loop(3,1),loop(3,2)
            do j=loop(2,1),loop(2,2)
                do i=loop(1,1),loop(1,2)
                    if(bc%bc_markers(i,j,k)==0)cycle
                    sumy=0.0_dp
                    do spec=1,species_number
                        Yloc(spec)=max(Y%pr(spec)%cells(i,j,k),0.0_dp)
                        sumy=sumy+Yloc(spec)
                    end do
                    if(sumy<=tiny(1.0_dp))cycle
                    Yloc=Yloc/sumy
                    M=this%thermo%thermo_ptr%mixture_molar_mass_from_mass_fractions(Yloc)
                    do spec=1,species_number
                        X(spec)=Yloc(spec)*M/this%thermo%thermo_ptr%molar_masses(spec)
                        Y%pr(spec)%cells(i,j,k)=Yloc(spec)
                    end do
                    h=this%thermo%thermo_ptr%mixture_enthalpy_molar(T%cells(i,j,k),X)
                    h0=this%thermo%thermo_ptr%mixture_enthalpy_molar(T_ref,X)
                    hs%cells(i,j,k)=(h-h0)/M
                    rho%cells(i,j,k)=this%thermodynamic_pressure*M/(r_gase_J*T%cells(i,j,k))
                end do
            end do
        end do
    end associate
    deallocate(X,Yloc)
end subroutine refresh_thermodynamic_ghosts

! Reconstruct the material characteristic family (T and Y_k) by CABARET
! extrapolation and source-shifted maximum-principle limiting.
subroutine reconstruct_contact_family_face_state(this)
    class(cabaret_low_mach_solver),intent(inout)::this

    integer::dimensions,species_number
    integer::i,j,k,dim,spec,bound_number
    integer,dimension(3,2)::inner_loop,outer_loop,loop
    real(dp)::T_face(2),T_new(2),T_limited(2),T_cell,dT_dt
    real(dp)::u_face
    real(dp),dimension(:,:),allocatable::Y_face,Y_new,Y_limited
    real(dp),dimension(:),allocatable::Y_cell,Y_state,rhoY_source,dY_dt

    dimensions=this%domain%get_domain_dimensions()
    species_number=this%chem%chem_ptr%species_number
    inner_loop=this%domain%get_local_inner_cells_bounds()
    outer_loop=this%domain%get_local_utter_cells_bounds()

    call this%exchange_face_pressure_density()
    call this%initialize_material_temperature_faces()

    associate(Y=>this%Y%v_ptr,T=>this%T%s_ptr,rho=>this%rho%s_ptr, &
            rhoY_source_old=>this%species_mass_source,species_face_old=>this%species_face_old, &
            T_f_new=>this%T_f_new%s_ptr,Y_f_new=>this%Y_f_new%v_ptr, &
            v_f_new=>this%v_f_new%v_ptr,bc=>this%boundary%bc_ptr)
        !$omp parallel default(shared) &
            !$omp& private(i,j,k,dim,spec,bound_number,loop,T_face,T_new,T_limited, &
            !$omp& T_cell,dT_dt,u_face,Y_face,Y_new,Y_limited,Y_cell,Y_state, &
            !$omp& rhoY_source,dY_dt)
        allocate(Y_face(species_number,2),Y_new(species_number,2), &
            Y_limited(species_number,2),Y_cell(species_number), &
            Y_state(species_number),rhoY_source(species_number), &
            dY_dt(species_number))

        do dim=1,dimensions
            loop(3,1)=outer_loop(3,1)*I_m(dim,3)+inner_loop(3,1)*(1-I_m(dim,3))
            loop(3,2)=outer_loop(3,2)*I_m(dim,3)+inner_loop(3,2)*(1-I_m(dim,3))
            loop(2,1)=outer_loop(2,1)*I_m(dim,2)+inner_loop(2,1)*(1-I_m(dim,2))
            loop(2,2)=outer_loop(2,2)*I_m(dim,2)+inner_loop(2,2)*(1-I_m(dim,2))
            loop(1,1)=outer_loop(1,1)*I_m(dim,1)+inner_loop(1,1)*(1-I_m(dim,1))
            loop(1,2)=outer_loop(1,2)*I_m(dim,1)+inner_loop(1,2)*(1-I_m(dim,1))

            !$omp do collapse(3) schedule(guided)
            do k=loop(3,1),loop(3,2)
                do j=loop(2,1),loop(2,2)
                    do i=loop(1,1),loop(1,2)
                        if (bc%bc_markers(i,j,k)/=0) cycle

                        do spec=1,species_number
                            Y_face(spec,1)=species_face_old(spec,dim,i,j,k)
                            Y_face(spec,2)=species_face_old(spec,dim,i+I_m(dim,1), &
                                j+I_m(dim,2),k+I_m(dim,3))
                            Y_cell(spec)=Y%pr(spec)%cells(i,j,k)
                            Y_state(spec)=Y_cell(spec)
                            rhoY_source(spec)=rhoY_source_old(spec,i,j,k)
                            Y_new(spec,1)=2.0_dp*Y_cell(spec)-Y_face(spec,2)
                            Y_new(spec,2)=2.0_dp*Y_cell(spec)-Y_face(spec,1)
                        end do

                        T_face(1)=T_f_new%cells(dim,i,j,k)
                        T_face(2)=T_f_new%cells(dim,i+I_m(dim,1), &
                            j+I_m(dim,2),k+I_m(dim,3))
                        T_cell=T%cells(i,j,k)
                        T_new(1)=2.0_dp*T_cell-T_face(2)
                        T_new(2)=2.0_dp*T_cell-T_face(1)

                        call this%material_characteristic_source_rates(rho%cells(i,j,k),T_cell,Y_state, &
                            this%sensible_enthalpy_source(i,j,k)+ &
                            merge(this%thermodynamic_pressure_rate_mid,0.0_dp,this%closed_domain), &
                            rhoY_source,dY_dt,dT_dt)

                        T_limited(1)=limit_quasi_invariant(T_new(1),T_face(1),T_cell,T_face(2), &
                            dT_dt*this%time_step, nearest(1.0_dp,1.0_dp), &
                            nearest(tables_temperature_ceiling,-1.0_dp))
                        T_limited(2)=limit_quasi_invariant(T_new(2),T_face(1),T_cell,T_face(2), &
                            dT_dt*this%time_step, nearest(1.0_dp,1.0_dp), &
                            nearest(tables_temperature_ceiling,-1.0_dp))
                        do spec=1,species_number
                            Y_limited(spec,1)=limit_quasi_invariant(Y_new(spec,1),Y_face(spec,1), &
                                Y_cell(spec),Y_face(spec,2),dY_dt(spec)*this%time_step, 0.0_dp, 1.0_dp)
                            Y_limited(spec,2)=limit_quasi_invariant(Y_new(spec,2),Y_face(spec,1), &
                                Y_cell(spec),Y_face(spec,2),dY_dt(spec)*this%time_step, 0.0_dp, 1.0_dp)
                        end do

#ifdef OMP
                        call omp_set_lock(lock(i,j,k))
                        call omp_set_lock(lock(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
#endif
                        if ((I_m(dim,1)*i+I_m(dim,2)*j+I_m(dim,3)*k)/=outer_loop(dim,1)) then
                            bound_number=bc%bc_markers(i-I_m(dim,1), &
                                j-I_m(dim,2),k-I_m(dim,3))
                            u_face=v_f_new%pr(dim)%cells(dim,i,j,k)
                            if (bound_number==0 .and. u_face<0.0_dp) then
                                do spec=1,species_number
                                    Y_f_new%pr(spec)%cells(dim,i,j,k)=Y_limited(spec,1)
                                end do
                                T_f_new%cells(dim,i,j,k)=T_limited(1)
                            end if
                        end if
                        if ((I_m(dim,1)*i+I_m(dim,2)*j+I_m(dim,3)*k)/=outer_loop(dim,2)) then
                            bound_number=bc%bc_markers(i+I_m(dim,1), &
                                j+I_m(dim,2),k+I_m(dim,3))
                            u_face=v_f_new%pr(dim)%cells(dim,i+I_m(dim,1), &
                                j+I_m(dim,2),k+I_m(dim,3))
                            if (bound_number==0 .and. u_face>=0.0_dp) then
                                do spec=1,species_number
                                    Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1), &
                                        j+I_m(dim,2),k+I_m(dim,3))=Y_limited(spec,2)
                                end do
                                T_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2), &
                                    k+I_m(dim,3))=T_limited(2)
                            end if
                        end if
#ifdef OMP
                        call omp_unset_lock(lock(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
                        call omp_unset_lock(lock(i,j,k))
#endif
                    end do
                end do
            end do
            !$omp end do
        end do

        deallocate(Y_face,Y_new,Y_limited,Y_cell,Y_state,rhoY_source,dY_dt)
        !$omp end parallel
    end associate
end subroutine reconstruct_contact_family_face_state

	subroutine finalize_gas_dynamics_step(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

        associate(	v   => this%v%v_ptr)
			call this%mpi_support%exchange_conservative_vector_field(v)
        end associate

		call this%apply_cell_boundary_conditions()
	end subroutine finalize_gas_dynamics_step

subroutine advance_particle_phases(this)
    class(cabaret_low_mach_solver), intent(inout) :: this
    integer :: phase

    do phase = 1, this%additional_particles_phases_number
        call this%particles_solver(phase)%advance(this%time_step, this%time)
    end do
end subroutine advance_particle_phases



	subroutine cache_face_state_for_next_step(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		integer :: dimensions, species_number
		integer :: i, j, k, dim, component, spec
		integer, dimension(3,2) :: face_loop

		dimensions = this%domain%get_domain_dimensions()
		species_number = this%chem%chem_ptr%species_number
		face_loop = this%domain%get_local_utter_faces_bounds()

		! CABARET face variables are independent degrees of freedom.  Preserve
		! the projected/reconstructed state as the old-time face state used by
		! the next predictor and characteristic extrapolation.
		!$omp parallel default(shared) private(i,j,k,dim,component,spec)
		do dim = 1, dimensions
			!$omp do collapse(3) schedule(static)
			do k = face_loop(3,1), face_loop(3,2)
			do j = face_loop(2,1), face_loop(2,2)
			do i = face_loop(1,1), face_loop(1,2)
				this%density_face_old(dim,i,j,k) = this%rho_f_new%s_ptr%cells(dim,i,j,k)
				do spec = 1, species_number
					this%species_face_old(spec,dim,i,j,k) = &
						this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k)
				end do
				do component = 1, dimensions
					this%velocity_face_old(component,dim,i,j,k) = &
						this%v_f_new%v_ptr%pr(component)%cells(dim,i,j,k)
				end do
			end do
			end do
			end do
			!$omp end do nowait
		end do
		!$omp end parallel
	end subroutine cache_face_state_for_next_step

	! State exchange and persistence helpers
	! The original implementation had several repeated OpenMP/MPI blocks inside
	! solve_problem().  The routines below collect the repeated state management
	! operations so that the numerical CABARET stages are easier to read and
	! modify.  They intentionally do not alter the mathematical algorithm except
	! for the explicitly marked consistency fixes in the face-species treatment.

	subroutine cache_cell_state(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		integer :: dim, spec, dimensions, species_number

		dimensions      = this%domain%get_domain_dimensions()
		species_number  = this%chem%chem_ptr%species_number

		this%rho_old = this%rho%s_ptr%cells
		this%h_s_old = this%h_s%s_ptr%cells

		do spec = 1, species_number
			this%Y_old(spec,:,:,:) = this%Y%v_ptr%pr(spec)%cells
		end do

		do dim = 1, dimensions
			this%v_old(dim,:,:,:) = this%v%v_ptr%pr(dim)%cells
		end do
	end subroutine cache_cell_state

	subroutine exchange_conservative_state(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		call this%mpi_support%exchange_conservative_scalar_field(this%p%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%rho%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%E_f%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%v_s%s_ptr)
		call this%mpi_support%exchange_conservative_scalar_field(this%h_s%s_ptr)

		call this%mpi_support%exchange_conservative_vector_field(this%Y%v_ptr)
		call this%mpi_support%exchange_conservative_vector_field(this%v%v_ptr)
	end subroutine exchange_conservative_state

	subroutine exchange_face_pressure_density(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
	end subroutine exchange_face_pressure_density

	subroutine exchange_face_primitive_state(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_scalar_field(this%rho_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%p_f_new%s_ptr)
		call this%mpi_support%exchange_flow_scalar_field(this%T_f_new%s_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_primitive_state

	subroutine exchange_face_thermodynamic_state(this)
		class(cabaret_low_mach_solver), intent(inout) :: this

		call this%mpi_support%exchange_flow_vector_field(this%v_f_new%v_ptr)
		call this%mpi_support%exchange_flow_vector_field(this%Y_f_new%v_ptr)
	end subroutine exchange_face_thermodynamic_state

subroutine initialize_material_temperature_faces(this)
    class(cabaret_low_mach_solver),intent(inout)::this
    integer::dimensions,dim,dim1,i,j,k
    integer,dimension(3,2)::face_loop,loop

    dimensions=this%domain%get_domain_dimensions()
    face_loop=this%domain%get_local_inner_faces_bounds()
    associate(density_face_old=>this%density_face_old,species_face_old=>this%species_face_old,T_f_new=>this%T_f_new%s_ptr)
    do dim=1,dimensions
        loop=face_loop
        do dim1=1,dimensions
            loop(dim1,2)=face_loop(dim1,2)-(1-I_m(dim1,dim))
        end do
        !$omp parallel do collapse(3) schedule(static) private(i,j,k)
        do k=loop(3,1),loop(3,2)
        do j=loop(2,1),loop(2,2)
        do i=loop(1,1),loop(1,2)
            T_f_new%cells(dim,i,j,k)=this%temperature_face_value( &
                density_face_old(dim,i,j,k),species_face_old(:,dim,i,j,k))
        end do
        end do
        end do
        !$omp end parallel do
    end do
    end associate
end subroutine initialize_material_temperature_faces

subroutine normalize_face_mass_fractions(this)
    class(cabaret_low_mach_solver), intent(inout) :: this

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

                    ! fall back to the old face state.  This avoids a zero molecular weight in
                    ! apply_state_equation_flow_variables().
                    if (spec_summ <= tiny(1.0_dp)) then
                        spec_summ = 0.0_dp
                        do spec = 1, species_number
                            this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max(this%species_face_old(spec,dim,i,j,k), 0.0_dp)
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
                                this%Y_f_new%v_ptr%pr(spec)%cells(dim,i,j,k) = max( &
                                    this%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)),0.0_dp)
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
	! Populate cell ghost states for wall, symmetry, inlet, outlet, and
	! backflow conditions; the projection boundary conditions are set separately.
	subroutine apply_cell_boundary_conditions(this)

		class(cabaret_low_mach_solver)		,intent(inout)		:: this

		integer					:: dimensions
		integer, dimension(3,2) :: cons_inner_loop
		character(len=20)		:: boundary_type_name
		real(dp)				:: farfield_density, farfield_pressure, farfield_temperature, wall_temperature, farfield_velocity
		real(dp)				:: spec_summ, mol_mix_from_farfield, normal_velocity

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,specie_number
		integer :: ghost_i, ghost_j, ghost_k

		dimensions			= this%domain%get_domain_dimensions()
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
										write(*,*) 'CABARET apply_cell_boundary_conditions: unsupported boundary type ', &
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

! Apply low-Mach advective and explicit parabolic stability restrictions;
! acoustic speed does not enter the projection-method time step.
subroutine calculate_time_step(this)
#ifdef mpi
    use MPI
#endif
    class(cabaret_low_mach_solver), intent(inout) :: this

    integer :: dimensions, species_number
    integer :: i, j, k, dim, spec
    integer, dimension(3,2) :: loop
    real(dp), dimension(3) :: dx
    real(dp) :: inverse_dx2_sum, advective_rate, diffusivity_max
    real(dp) :: dt_advective_local, dt_parabolic_local, dt_local, dt_global
    real(dp) :: thermal_diffusivity
#ifdef mpi
    integer :: error, mpi_communicator
#endif

    dimensions = this%domain%get_domain_dimensions()
    species_number = this%chem%chem_ptr%species_number
    loop = this%domain%get_local_inner_cells_bounds()
    dx = this%mesh%mesh_ptr%get_cell_edges_length()

    inverse_dx2_sum = 0.0_dp
    do dim = 1, dimensions
        if (dx(dim) <= 0.0_dp) error stop 'CABARET projection: invalid mesh spacing'
        inverse_dx2_sum = inverse_dx2_sum+1.0_dp/dx(dim)**2
    end do

    dt_advective_local = huge(1.0_dp)
    dt_parabolic_local = huge(1.0_dp)

    associate(v => this%v%v_ptr, div_v => this%div_v%s_ptr, rho => this%rho%s_ptr, &
            cp => this%mixture_cp%s_ptr, M => this%mix_mol_mass%s_ptr, &
            bc => this%boundary%bc_ptr)

        do k = loop(3,1), loop(3,2)
            do j = loop(2,1), loop(2,2)
                do i = loop(1,1), loop(1,2)
                    if (bc%bc_markers(i,j,k) /= 0) cycle
                    if (rho%cells(i,j,k) <= 0.0_dp) &
                        error stop 'CABARET projection: non-positive density in time-step estimate'

                    advective_rate = abs(div_v%cells(i,j,k))
                    do dim = 1, dimensions
                        advective_rate = advective_rate+abs(v%pr(dim)%cells(i,j,k))/dx(dim)
                    end do
                    if (advective_rate > 0.0_dp) &
                        dt_advective_local = min(dt_advective_local,1.0_dp/advective_rate)

                    diffusivity_max = 0.0_dp
                    if (this%viscosity_flag) &
                        diffusivity_max = max(diffusivity_max,this%nu%s_ptr%cells(i,j,k)/rho%cells(i,j,k))
                    if (this%heat_trans_flag) then
                        if (cp%cells(i,j,k) <= 0.0_dp .or. M%cells(i,j,k) <= 0.0_dp) &
                            error stop 'CABARET projection: invalid heat capacity in time-step estimate'
                        thermal_diffusivity = this%kappa%s_ptr%cells(i,j,k)/ &
                            (rho%cells(i,j,k)*cp%cells(i,j,k)/M%cells(i,j,k))
                        diffusivity_max = max(diffusivity_max,thermal_diffusivity)
                    end if
                    if (this%diffusion_flag) then
                        do spec = 1, species_number
                            diffusivity_max = max(diffusivity_max,this%D%v_ptr%pr(spec)%cells(i,j,k))
                        end do
                    end if
                    if (diffusivity_max > 0.0_dp) &
                        dt_parabolic_local = min(dt_parabolic_local,0.25_dp/(diffusivity_max*inverse_dx2_sum))
                end do
            end do
        end do
    end associate

    dt_local = min(this%courant_fraction*dt_advective_local,dt_parabolic_local)
#ifdef mpi
    mpi_communicator = this%domain%get_mpi_communicator()
    call mpi_allreduce(dt_local,dt_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,mpi_communicator,error)
#else
    dt_global = dt_local
#endif
    if (dt_global < huge(1.0_dp)) then
        if (this%time_step > 0.0_dp) &
            dt_global = min(dt_global,cabaret_time_step_growth_limit*this%time_step)
        this%time_step = dt_global
    end if
end subroutine calculate_time_step

	subroutine set_CFL_coefficient(this,coefficient)
		class(cabaret_low_mach_solver)	,intent(inout)	:: this
		real(dp)				,intent(in)		:: coefficient

		this%courant_fraction = coefficient

	end subroutine

	pure function get_time_step(this)
		real(dp)						:: get_time_step
		class(cabaret_low_mach_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dp)						:: get_time
		class(cabaret_low_mach_solver)	,intent(in)		:: this

		get_time = this%time
	end function

end module
